from Bio import SeqIO
import pandas as pd
out_records = []
count = 1
mapping = []
with open(snakemake.input[0], 'r') as handle:
  for record in SeqIO.parse(handle, 'fasta'):
    if len(record.seq) >= snakemake.config["min_contig_size"]:
      new_id = f'{snakemake.param.prefix}_{count:07d}'
      record.id = f'{snakemake.param.prefix}_{count:07d}'
      mapping.append((record.id, new_id))
      record.id = new_id
      record.description =''
      out_records.append(record)

SeqIO.write(out_records, snakemake.output[0], 'fasta')
pd.DataFrame(mapping, columns=['old_id', 'new_id']).to_csv(snakemake.output[1], sep='\t', index=False)

      
