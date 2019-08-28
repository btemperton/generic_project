guppy_version=['guppy_v322']

rule all:
    input:
        expand('/dbs/wtg/{sample}.{guppy_version}.fast5.tar.gz',
                sample=samples, guppy_version=guppy_version),
        expand('/dbs/wtg/{sample}.{guppy_version}.fast5.tar.gz.md5',
                sample=samples, guppy_version=guppy_version),
        expand('/dbs/wtg/{sample}.{guppy_version}.pass.fq.gz',
                sample=samples, guppy_version=guppy_version),
        expand('/dbs/wtg/{sample}.{guppy_version}.sequencing_summary.txt.gz',sample=samples, guppy_version=guppy_version)

rule basecall:
    input: "/reads/projects/TempertonLab/wtg/samples/{sample}/reads/fast5"
    output:
        temp(directory('/dbs/wtg/{sample}.{guppy_version}'))
    params:
        device='cuda:0',
        config_file='dna_r9.4.1_450bps_hac.cfg'
    log:  'logs/{sample}.{guppy_version}.basecall.log'
    benchmark:  'benchmark/{sample}.{guppy_version}.tsv'
    shell:
        """
            guppy_basecaller -i {input} \
            -s {output} \
            --device {params.device} \
            -c {params.config_file} \
            --trim_strategy dna \
            --records_per_fastq 0 \
            --qscore_filtering \
            --min_qscore 7 \
            --fast5_out 2>&1 | tee {log};
        """

rule extract_pass:
    input: rules.basecall.output
    output: '/dbs/wtg/{sample}.{guppy_version}.pass.fq.gz'
    shell:
        """
            cat {input}/pass/*.fastq > /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.pass.fq;
            pigz /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.pass.fq
        """

rule extract_fast5:
    input: rules.basecall.output
    output: "/dbs/wtg/{sample}.{guppy_version}.fast5.tar.gz"
    shell:
        """
        find {input} -name "*.fast5" | tar cf /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.fast5.tar -T -;
        pigz /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.fast5.tar;
        """
rule create_md5:
    input: rules.extract_fast5.output
    output: "/dbs/wtg/{sample}.{guppy_version}.fast5.tar.gz.md5"
    shell:
        """
            md5sum {input} > {output}
        """

rule get_sequencing_summary:
    input: rules.basecall.output
    output: "/dbs/wtg/{sample}.{guppy_version}.sequencing_summary.txt.gz"
    shell:
        """
            cp {input}/sequencing_summary.txt /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.sequencing_summary.txt;
            pigz /dbs/wtg/{wildcards.sample}.{wildcards.guppy_version}.sequencing_summary.txt;
        """
rule long_qc:
    input: rules.extract_pass.output
    output:
        "/dbs/wtg/{sample}.{guppy_version}.porechop.1k5.fq.gz"
    threads: 16
    conda: "../envs/long_qc.yml"
    params:
        filter_quality=7,
        min_length=1500
    shell:
        """
        porechop -i {input} \
        --format fastq.gz \
        -t {threads} | \
        NanoFilt -q {params.filter_quality} \
        -l {params.min_length} \
        --headcrop 50 \
        --tailcrop 50 \
        --readtype 1D | gzip > {output}
        """
