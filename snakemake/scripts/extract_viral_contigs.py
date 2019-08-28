from Bio import SeqIO
import pandas as pd
import re
import glob

rgx = re.compile(r'(VIRSorter.*)(gene_\d+)_(gene_\d+)-\d+-\d+-cat_\d+')

def read_virsorter_data(output_dir, min_viral_length=5000):
    df = pd.read_csv(f'{output_dir}/VIRSorter_global-phage-signal.csv',
        header=None,
        usecols=[0,1,2,3,4],
        names=['contig_id', 'nb_genes_contigs', 'fragment', 'nb_genes', 'category'],
        comment='#')

    df['is_circular'] = df.contig_id.str.contains('circular')
    df['fragment'] = df.fragment.str.replace('-circular-gene', '-circular_gene')
    repl = lambda m: m.group(1)
    df['original_contig'] = df.contig_id.str.replace('VIRSorter_', '')
    df['original_contig'] = df.original_contig.str.replace('-circular', '')
    reads = {}
    lengths = []
    for f in glob.glob(f'{output_dir}/Predicted_viral_sequences/*.fasta'):
        with open(f, 'r') as handle:
            for r in SeqIO.parse(handle, 'fasta'):
                m = rgx.search(r.id)
                if m:
                    r.id = f'{m.group(1)}{m.group(2)}-{m.group(3)}'
                else:
                    r.id = r.id.split('-cat_')[0]
                reads[r.id] = r
                lengths.append((r.id, len(r.seq)))

    length_df = pd.DataFrame(lengths, columns=['fragment', 'contig_length'])

    new_df = df.merge(length_df, how='left', on='fragment')

    vs_viral_df = new_df[(new_df.is_circular) | (new_df.contig_length >= min_viral_length)]

    return vs_viral_df, reads

df, reads = read_virsorter_data('.')
