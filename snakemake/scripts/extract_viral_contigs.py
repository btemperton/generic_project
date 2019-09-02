from Bio import SeqIO
import pandas as pd
import re
import glob
import numpy as np

rgx = re.compile(r'(VIRSorter.*)(gene_\d+)_(gene_\d+)-\d+-\d+-cat_\d+')

def read_virsorter_data(output_dir):
    df = pd.read_csv(f'{output_dir}/VIRSorter_global-phage-signal.csv',
        header=None,
        usecols=[0,1,2,3,4],
        names=['contig_id', 'nb_genes_contigs', 'fragment', 'nb_genes', 'category'],
        comment='#')

    df.loc[df.fragment.str.contains(r'-gene_\d+-gene_\d+'), 'category'] += 3

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

    length_df = pd.DataFrame(lengths, columns=['fragment', 'vs_contig_length'])

    new_df = df.merge(length_df, how='left', on='fragment')

    return new_df, reads

def read_virfinder_output(vf_output_file):
    if vf_output_file:
        df = pd.read_csv(vf_output_file, sep='\t',
            names=['original_contig','vf_contig_length', 'vf_score', 'vf_pvalue'], skiprows=1)
    else:
        df = None
    return df


def get_contig_lengths(contig_file):
    original_contigs = {}
    lengths = []
    with open(contig_file, 'r') as handle:
        for r in SeqIO.parse(handle, 'fasta'):
            original_contigs[r.id] = r
            lengths.append((r.id, len(r.seq)))
    return pd.DataFrame(lengths, columns=['original_contig', 'length']), original_contigs


def extract_viral_contigs(contig_file, vs_output, vf_output, combined_output, fasta_output, min_contig_length=5000):
    lengths_df, original_contigs = get_contig_lengths(contig_file)
    vs_df, vs_contigs =  read_virsorter_data(vs_output)
    vf_df = read_virfinder_output(vf_output)

    new_df = lengths_df.merge(vs_df, how='left', on='original_contig')
    if vf_df is not None:
        new_df = new_df.merge(vf_df, how='left', on='original_contig')

    conditions = [((new_df.category <3) & ((new_df.vs_contig_length >= min_contig_length) | (new_df.is_circular))),
                    ((new_df.vs_contig_length >= min_contig_length) | (new_df.is_circular))]
    choices = ['Definite', 'Putative']

    new_df['is_vs_viral'] = np.select(conditions, choices, default='Not')

    if vf_df is not None:
        conditions = [((new_df.vf_score >= 0.9) & (new_df.vf_pvalue <0.05)),
                        ((new_df.vf_score >= 0.7) & (new_df.vf_pvalue <0.05))]
        new_df['is_vf_viral'] = np.select(conditions, choices, default='Not')

    conditions = [(new_df.is_vs_viral == 'Definite'),
                  (vf_df is not None and new_df.is_vf_viral == 'Definite'),
                  ((new_df.is_vs_viral == 'Putative') & (vf_df is not None and new_df.is_vf_viral == 'Putative'))]

    choices = ['vs_viral', 'vf_viral', 'combined_viral']

    new_df['viral_status'] = np.select(conditions, choices, default='Not')
    new_df.to_csv(combined_output, sep='\t', index=False)

    contigs_out = []

    for i in new_df[new_df.viral_status == 'vs_viral']['contig_id'].values:
        contigs_out.append(vs_contigs[i])
    if vf_df is not None:
        for i in new_df[new_df.viral_status == 'vf_viral']['original_contig'].values:
            contigs_out.append(original_contigs[i])
        for i in new_df[new_df.viral_status == 'combined_viral']['original_contig'].values:
            contigs_out.append(original_contigs[i])

    SeqIO.write(contigs_out, fasta_output, 'fasta')

extract_viral_contigs(snakemake.input.contigs,
                        snakemake.input.virsorter_out,
                        snakemake.input.virfinder_out,
                        snakemake.output.viral_mapping,
                        snakemake.output.viral_contigs,
                        snakemake.params.min_non_circular_viral_length)
