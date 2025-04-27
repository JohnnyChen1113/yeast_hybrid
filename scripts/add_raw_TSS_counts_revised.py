import click
import pandas as pd
import numpy as np
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@click.command()
@click.option('-c', '--clusters-file', required=True, type=click.Path(exists=True), help='Path to assignedClusters file')
@click.option('-t', '--tss-file', required=True, type=click.Path(exists=True), help='Path to raw TSS file')
@click.option('--merge', is_flag=True, help='Merge TSS counts from all samples into a single column')
@click.option('-o', '--output-file', required=True, type=click.Path(), help='Path to output file')
@click.option('--chunk-size', default=100000, type=int, help='Chunk size for reading TSS file')
def add_tss_counts_optimized(clusters_file, tss_file, merge, output_file, chunk_size):
    """
    Efficiently add TSS counts from raw TSS file to assignedClusters file for large-scale data.
    
    This script processes large assignedClusters and raw TSS files, calculating TSS counts within
    cluster regions (based on chr, strand, and pos within start-end range) and appends them as new columns.
    Optimized for performance using vectorized operations and chunked file reading.
    """
    logging.info("Starting processing...")
    
    # Read assignedClusters file
    logging.info(f"Reading clusters file: {clusters_file}")
    clusters_df = pd.read_csv(clusters_file, sep='\t')
    
    # Ensure clusters_df has a unique index
    clusters_df = clusters_df.reset_index(drop=True)
    
    # Read raw TSS file to get sample columns
    logging.info(f"Reading TSS file header: {tss_file}")
    tss_header = pd.read_csv(tss_file, sep='\t', nrows=0)
    sample_columns = tss_header.columns[3:]
    logging.info(f"Detected sample columns: {sample_columns}")
    
    # Initialize count columns in clusters_df
    logging.info("Initializing count columns in clusters DataFrame")
    if merge:
        clusters_df['merged_TSS_counts'] = 0
    else:
        for col in sample_columns:
            clusters_df[f'{col}_TSS_counts'] = 0
    
    # Process TSS file in chunks
    logging.info(f"Processing TSS file in chunks of size {chunk_size}")
    for i, chunk in enumerate(pd.read_csv(tss_file, sep='\t', chunksize=chunk_size)):
        logging.info(f"Processing chunk {i + 1}")
        
        # Add original cluster index to clusters_df for merging
        clusters_df['cluster_index'] = clusters_df.index
        
        # Merge clusters_df with chunk on chr and strand
        merged = pd.merge(
            clusters_df[['cluster_index', 'chr', 'start', 'end', 'strand']],
            chunk,
            on=['chr', 'strand'],
            how='inner'
        )
        
        # Filter rows where pos is within [start, end]
        valid = (merged['pos'] >= merged['start']) & (merged['pos'] <= merged['end'])
        merged = merged[valid]
        
        if not merged.empty:
            # Group by original cluster_index
            if merge:
                # Sum all sample columns
                merged['merged_TSS_counts'] = merged[sample_columns].sum(axis=1)
                counts = merged.groupby('cluster_index')['merged_TSS_counts'].sum()
                clusters_df.loc[counts.index, 'merged_TSS_counts'] += counts
            else:
                # Sum each sample column separately
                for col in sample_columns:
                    counts = merged.groupby('cluster_index')[col].sum()
                    clusters_df.loc[counts.index, f'{col}_TSS_counts'] += counts
        else:
            logging.info("No valid TSS positions found in this chunk")
        
        # Clean up
        clusters_df = clusters_df.drop(columns=['cluster_index'])
    
    # Ensure count columns are integers
    count_cols = ['merged_TSS_counts'] if merge else [f'{col}_TSS_counts' for col in sample_columns]
    clusters_df[count_cols] = clusters_df[count_cols].astype(int)
    
    # Write to output file
    logging.info(f"Writing output to: {output_file}")
    clusters_df.to_csv(output_file, sep='\t', index=False)
    logging.info("Processing complete")

if __name__ == '__main__':
    add_tss_counts_optimized()
