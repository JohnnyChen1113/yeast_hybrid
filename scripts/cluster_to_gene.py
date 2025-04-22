import click
import pandas as pd

__version__ = "0.2.1"

@click.command()
@click.option('-i', '--input_file', type=click.Path(exists=True), required=True, help='Path to the input file')
@click.option('-o', '--output_file', type=click.Path(), help='Path to the output file, defaults to input file prefix with .gene.txt')
@click.option('-d', '--dominate', is_flag=True, help='Extract rows with the maximum value in the tags column as dominate')
def cluster_to_gene(input_file, output_file, dominate):
    """
    Converts a cluster-centric file to a gene-centric file, placing cluster_number in the second column.
    Optionally extracts rows with the maximum tags value when --dominate is specified.
    """
    try:
        df = pd.read_csv(input_file, sep='\s+')
    except FileNotFoundError:
        click.echo(f"Error: File not found at {input_file}")
        return
    except Exception as e:
        click.echo(f"Error reading file: {e}")
        return

    if 'gene' not in df.columns:
        click.echo("Error: Input file missing 'gene' column.")
        return

    if dominate:
        if 'tags' not in df.columns:
            click.echo("Error: Input file missing 'tags' column required for --dominate option.")
            return
        # Extract rows with the maximum tags value per gene
        df = df.loc[df.groupby('gene')['tags'].idxmax()]

    gene_grouped = df.groupby('gene').agg(lambda x: ','.join(map(str, x)))
    gene_grouped['cluster_number'] = df.groupby('gene').size()
    gene_grouped.reset_index(inplace=True)

    # Adjust the order of columns
    columns = gene_grouped.columns.tolist()
    columns.insert(1, columns.pop(columns.index('cluster_number')))
    gene_grouped = gene_grouped[columns]

    # Replace nan with NA in the output
    gene_grouped = gene_grouped.fillna('NA')

    if output_file:
        gene_grouped.to_csv(output_file, sep='\t', index=False)
        click.echo(f"Results saved to: {output_file}")
    else:
        output_file_prefix = input_file.rsplit('.', 1)[0] if '.' in input_file else input_file
        default_output_file = f"{output_file_prefix}.gene.txt"
        gene_grouped.to_csv(default_output_file, sep='\t', index=False)
        click.echo(f"Results saved to default file: {default_output_file}")

if __name__ == '__main__':
    cluster_to_gene()
