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

rule extract_viral_contigs:
    input:
        contigs=rules.rename_contigs.output.contigs,
        virsorter_out=rules.virsorter.output,
        virfinder_out=rules.DeepVirFinder.output
    output:
        viral_contigs="samples/{sample}/assembly/viral_contigs.fa",
        viral_mapping="samples/{sample}/assembly/viral_id.tsv"
    params:
        min_non_circular_viral_length=5000
    log:
        "logs/{sample}.extract_viral_contigs.log"
    benchmark:
        "benchmarks/{sample}.extract_viral_contigs.tsv"
    script:
        "../scripts/extract_viral_contigs.py"
