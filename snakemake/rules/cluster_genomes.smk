rule cluster_genomes:
  input:
    "something goes here"
  output:
    "something goes here:
  conda:
        "../envs/cluster_genomes.yml"
    params:
        coverage=80,
        identity=95
    log:
        "logs/{sample}.cluster_genomes.log"
    benchmark:
        "benchmarks/{sample}.cluster_genomes.tsv"
    script:
      "../scripts/cluster_genomes.py"
