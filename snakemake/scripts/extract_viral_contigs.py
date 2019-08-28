from Bio import SeqIO
import pandas as pd
import re
import glob


def read_virsorter_data(output_dir, min_viral_length):
    df = pd.read_csv(f'{output_dir}/VIRSorter_global-phage-signal.csv',
        header=None,
        usecols=[0,1,2,3,4],
        names=['contig_id', 'nb_genes_contigs', 'fragment', 'nb_genes', 'category'],
        comment='#')

    df['is_circular'] = df.contig_id.str.contains('circular')

    reads = {}
    lengths = []
    for f in glob.glob(f'{output_dir}/Predicted_viral_sequences/*.fasta'):
        with open(f, 'r') as handle:
            for r in SeqIO.parse(handle, 'fasta'):
                r.id = r.id.split('-cat_')[0]
                reads[r.id] = r
                lengths.append((r.id, len(r.seq)))

    length_df = pd.DataFrame(lengths, columns=['contig_id', 'contig_length'])

    new_df = df.merge(length_df, how='left', on='contig_id')

    vs_viral_df = new_df[(new_df.is_circular) | (new_df.contig_length >= min_viral_length)]

    ##FINDING LENGTHS OF prophage isn't working.

    return vs_viral_df, reads


def extract_viral_contigs(visorter_output_dir, input_contigs):
