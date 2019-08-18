from Bio import SeqIO
import pandas as pd

def rename_fastas(fasta_in, prefix, fasta_out, min_contig_size, mapping_file):
  out_records = []
  count = 1
  mapping = []
  with open(fasta_in, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
      if len(record.seq) >= min_contig_size:
        new_id = f'{prefix}_{count:07d}'
        mapping.append((record.id, new_id))
        record.id = new_id
        record.description =''
        out_records.append(record)
        count +=1

  SeqIO.write(out_records, fasta_out, 'fasta')
  pd.DataFrame(mapping, columns=['old_id', 'new_id']).to_csv(mapping_file, sep='\t', index=False)

rename_fastas(snakemake.input[0], 
              snakemake.param.prefix, 
              snakemake.output[0], 
              int(snakemake.config["min_contig_size"]), 
              snakemake.output[1])

