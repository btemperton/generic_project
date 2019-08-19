import pandas as pd
from common import BufferedOutputHandler
from Bio import SeqIO

def rename_fastas(fasta_in, prefix, fasta_out, min_contig_size, mapping_file, log_file):
  out_handler = BufferedOutputHandler(fasta_out, log_file)
  count = 1
  mapping = []
  with open(fasta_in, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
      if len(record.seq) >= min_contig_size:
        new_id = f'{prefix}_{count:07d}'
        mapping.append((record.id, new_id, len(record.seq)))
        record.id = new_id
        record.description =''
        record.name = ''
        out_handler.add_record(record)
        count +=1
  out_handler.close_out()
  pd.DataFrame(mapping, columns=['old_id', 'new_id', 'length']).to_csv(mapping_file, sep='\t', index=False)

rename_fastas(snakemake.input[0],
              snakemake.params.prefix,
              snakemake.output[0],
              int(snakemake.config["min_contig_size"]),
              snakemake.output[1],
              str(snakemake.log))


