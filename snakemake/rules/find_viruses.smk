rule virsorter:
  input:
    virsorter_home=rules.initialise_virsorter.output.home,
    data_dir=rules.initialse_virsorter.output.data,
    contigs=rules.rename_contigs.output.contigs
  output:
    dir("samples/{sample}/assembly/virsorter")
  conda: "../envs/virsorter.yml"
  threads: 20
  shell:
    """
    {input.virsorter_home}/wrapper_phage_contigs_sorter_iPlant.pl \
    --fna {input.contigs} --db 2 --ncpu {threads} --virome --data-dir {input.data_dir} --diamond
    """
    
rule DeepVirFinder:
  input:
    contigs=rules.rename_contigs.output.contigs,
    dvf=rules.initialise_DeepVirFinder.cmd,
    dvf_home=rules.initialise_DeepVirFinder.directory
  output:
    dir("samples/{sample}/assembly/dvf")
  conda: "../envs/DeepVirFinder.yml"
  threads: 20
  shell:
    """
    python {input.dvf_home}/dvf.py -i {input.contigs} \
    -m {input.dvf_home}/models \
    -o {output} \
    -c {threads}
    """
