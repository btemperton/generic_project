import os
import pandas as pd
configfile: "config.yaml"

samples = pd.read_csv(config['sample_file']).set_index("sample", drop=False)
def get_fwd_url(wildcard):
    return samples.loc[wildcard, 'fwd'].values[0]

def get_rev_url(wildcard):
    return samples.loc[wildcard, 'rev'].values[0]

localrules: all,get_fwd_reads,get_rev_reads

rule all:
    input:
        expand("samples/{sample}/fwd.hq.gz",
                sample=samples.index),
        expand("samples/{sample}/rev.hq.gz",
                sample=samples.index)

rule get_fwd_reads:
    output:
        temp("samples/{sample}/fwd.gz")
    threads: 1
    params:
        fwd_url=get_fwd_url,
    log:
        "logs/{sample}.get_fwd_reads.log"
    benchmark:
        "benchmarks/{sample}.get_fwd_reads.tsv"
    shell:
        """
        wget -O {output} {params.fwd_url}
        """

rule get_rev_reads:
    output:
        temp("samples/{sample}/rev.gz")
    threads: 1
    params:
        rev_url=get_rev_url,
    log:
        "logs/{sample}.get_rev_reads.log"
    benchmark:
        "benchmarks/{sample}.get_rev_reads.tsv"
    shell:
        """
        wget -O {output} {params.rev_url}
        """
        
rule run_bbmap_clumpify:
    input:
        raw_fwd=rules.get_fwd_reads.output,
        raw_rev=rules.get_rev_reads.output
    output:
        temp("{sample}.clumped.fq.gz")
    threads: 16
    conda:
        "envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_clumpify.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_clumpify.tsv"
    group: "bbtools"
    shell:
        """
            clumpify.sh -Xmx104g -eoom -da in1={input.raw_fwd} in2={input.raw_rev} out={output} dedupe optical
        """

rule run_bbmap_filter_by_tile:
    input:
        rules.run_bbmap_clumpify.output
    output:
        temp("{sample}.filtered_by_tile.fq.gz")
    threads: 16
    conda:
        "envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_filter_by_tile.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_filter_by_tile.tsv"
    group: "bbtools"
    shell:
        """
            filterbytile.sh -eoom -Xmx104g -da in={input} out={output}
        """

rule run_bbmap_bbduk_remove_adapters:
    input:
        rules.run_bbmap_filter_by_tile.output
    output:
        temp("{sample}.trimmed.fq.gz")
    threads: 16
    conda:
        "envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_bbduk_remove_adapters.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_bbduk_remove_adapters.tsv"
    group: "bbtools"
    shell:
        """
            bbduk.sh -Xmx104g -eoom -da in={input} out={output} ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=50 ref=adapters ftm=5 ordered
        """

rule run_bbmap_bbduk_remove_artefacts:
    input:
        rules.run_bbmap_bbduk_remove_adapters.output
    output:
        temp("{sample}.filtered.fq.gz")
    threads: 16
    conda:
        "envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_bbduk_remove_artefacts.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_bbduk_remove_artefacts.tsv"
    group: "bbtools"
    shell:
        """
            bbduk.sh -Xmx104g -eoom -da in={input} out={output} k=31 ref=artifacts,phix ordered cardinality
        """

rule split_reads:
    input:
        rules.run_bbmap_bbduk_remove_artefacts.output
    output:
        "samples/{sample}/fwd.hq.gz",
        "samples/{sample}/rev.hq.gz"
    threads: 16
    conda:
        "envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.split_reads.log"
    benchmark:
        "benchmarks/{sample}.split_reads.tsv"
    params:
        sample_seed=42
    group: "bbtools"
    shell:
        """
            reformat.sh -Xmx104g -eoom -da in={input} out1={output[0]} out2={output[1]} sampleseed={params.sample_seed}
        """
