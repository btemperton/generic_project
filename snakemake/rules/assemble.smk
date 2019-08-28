rule spades_meta:
    input:
        fwd='samples/{sample}/reads/fwd.ec.hq.fq.gz',
        rev='samples/{sample}/reads/rev.ec.hq.fq.gz'
    output:
        contigs='samples/{sample}/assembly/spades.contigs.fa',
        fastg='samples/{sample}/assembly/spades.contigs.fastg',
    conda:
        "../envs/spades_3.13.1.yml"
    threads: 32
    params:
        memory=1024
    log:
        "logs/{sample}.spades_meta.log"
    benchmark:
        "benchmarks/{sample}.spades_meta.tsv"
    shell:
        """
        spades.py -o {wildcards.sample}.spades --meta \
        -1 {input.fwd} \
        -2 {input.rev} \
        --threads {threads} \
        --only-assembler \
        --memory {params.memory} 2>&1 | tee -a {log};

        mv {wildcards.sample}.spades/contigs.fasta {output.contigs};
        mv {wildcards.sample}.spades/assembly_graph.fastg {output.fastg};
        rm -rf {wildcards.sample}.spades;
        """

rule rename_contigs:
    input: rules.spades_meta.output.contigs
    output:
        contigs='samples/{sample}/assembly/spades.contigs.renamed.fa',
        mapping='samples/{sample}/assembly/mapping.tsv'
    params:
        prefix="{sample}"
    log:
        "logs/{sample}.rename_contigs.log"
    benchmark:
        "benchmarks/{sample}.rename_contigs.tsv"
    script:
        "../scripts/rename_fastas.py"
        
 rule bowtie2_build:
    input: rules.rename_contigs.output.contigs
    output: temp(dir('samples/{sample}/assembly/bowtie_idx'))
    threads: 20
    log:
        "logs/{sample}.bowtie2-build.log"
    benchmark:
        "benchmarks/{sample}.bowtie2-build.tsv"
    shell:
     """
     bowtie2-build --threads {threads} {input} {output}/bt
     """

rule bowtie2_map:
    input:
        bt_index=rules.bowtie2_build.output,
        fwd='samples/{sample}/reads/fwd.ec.hq.fq.gz',
        rev='samples/{sample}/reads/rev.ec.hq.fq.gz'
    output:
        'samples/{sample}/assembly/mapped.sorted.bam'
    threads: 20
    log:
        "logs/{sample}.bowtie2-build.log"
    benchmark:
        "benchmarks/{sample}.bowtie2-build.tsv"
    shell:
     """
     bowtie2 --threads {threads} -x {input.bt_index}/bt {input.fwd} {input.rev} | \
     samtools view -@ {threads} -F 4 -buS | \
     samtools sort -@ {threads} -o {output};
     samtools index {output}
     """
