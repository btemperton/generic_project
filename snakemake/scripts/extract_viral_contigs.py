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

    new_df['vs_viral'] = new_df.is_circular or
                        new_df.contig_length >= min_viral_length

    return new_df, reads

def read_virfinder_output(vf_output_file, min_score=0.7, max_p=0.05):
    if vf_output_file:
        df = pd.read_csv(vf_output, sep='\t',
            names=['original_contig','vf_contig_length', 'vf_score', 'vf_pvalue'], skiprows=1)
            df['vf_viral'] = (df.vf_score >= min_score) & (df.vf_pvalue <= max_p)
    else:
        df = None
    return df


def get_contig_lengths(contig_file):
    lengths = [(r.id, len(r.seq)) for r in SeqIO.parse(contig_file, 'fasta')]
    return pd.DataFrame(lengths, columns=['original_contig', 'length'])


def extract_viral_contigs(contig_file, vs_output, vf_output):
    lengths_df = get_contig_lengths(contig_file)
    vs_df, vs_contigs =  read_virsorter_data(vs_output)
    vf_df = read_virfinder_output(vf_output)
