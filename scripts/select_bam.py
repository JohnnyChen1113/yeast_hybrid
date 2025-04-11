import click
import pysam
import os

@click.command()
@click.option('-i', '--input', required=True, type=click.Path(exists=True), help='Input BAM file')
@click.option('-o', '--output', required=True, type=click.Path(), help='Output BAM file')
@click.option('-s', '--suffix', required=True, type=str, help='Two chromosome suffixes separated by comma (e.g., "_A,_B")')
def filter_bam(input, output, suffix):
    """
    Filter BAM file based on specific criteria:
    1. Filter out reads with flag=4
    2. Keep only reads with NH:i:1 or NH:i:2
    3. For NH:i:2, check if both suffixes have one record each, modify quality if true
    """
    # Check if index file exists, create it if not
    index_file = input + '.bai'
    if not os.path.exists(index_file):
        click.echo(f"Index file {index_file} not found. Generating index...")
        pysam.index(input)
        click.echo("Index generated successfully.")

    # Split suffixes
    suffix1, suffix2 = suffix.split(',')
    
    # Open input and output BAM files
    with pysam.AlignmentFile(input, 'rb') as in_bam, \
         pysam.AlignmentFile(output, 'wb', template=in_bam) as out_bam:
        
        # Dictionary to store NH:i:2 reads temporarily
        nh2_reads = {}
        
        # First pass: process reads
        for read in in_bam.fetch():
            # Filter out unmapped reads (flag=4)
            if read.flag == 4:
                continue
                
            # Get NH tag
            try:
                nh_value = read.get_tag('NH')
            except KeyError:
                continue
                
            # Process NH:i:1 reads - write directly
            if nh_value == 1:
                out_bam.write(read)
                continue
                
            # Process NH:i:2 reads - store temporarily
            if nh_value == 2:
                read_name = read.query_name
                if read_name not in nh2_reads:
                    nh2_reads[read_name] = []
                nh2_reads[read_name].append(read)
        
        # Second pass: process NH:i:2 reads
        for read_name, reads in nh2_reads.items():
            if len(reads) != 2:  # Must have exactly 2 reads
                continue
                
            # Check if reads match both suffixes
            chroms = [read.reference_name for read in reads]
            has_suffix1 = any(chrom.endswith(suffix1) for chrom in chroms)
            has_suffix2 = any(chrom.endswith(suffix2) for chrom in chroms)
            
            if has_suffix1 and has_suffix2:
                # Sort reads by flag value
                reads.sort(key=lambda x: x.flag)
                # Modify quality of read with smaller flag to 59
                reads[0].mapping_quality = 59
                # Write both reads
                for read in reads:
                    out_bam.write(read)

if __name__ == '__main__':
    filter_bam()
