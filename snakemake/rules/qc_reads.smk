def get_fwd_url(wildcard):
    return samples.loc[wildcard, 'fwd'].values[0]

def get_rev_url(wildcard):
    return samples.loc[wildcard, 'rev'].values[0]

rule get_reads:
    output:
        fwd=temp("samples/{sample}/fwd.gz"),
        rev=temp("samples/{sample}/rev.gz")
    threads: 1
    params:
        fwd_url=get_fwd_url,
        rev_url=get_rev_url
    log:
        "logs/{sample}.get_reads.log"
    benchmark:
        "benchmarks/{sample}.get_reads.tsv"
    shell:
        """
        wget -O {output.fwd} {params.fwd_url};
        wget -O {output.rev} {params.rev_url};
        """

rule run_bbmap_clumpify:
    input:
        raw_fwd=rules.get_reads.output.fwd,
        raw_rev=rules.get_reads.output.rev
    output:
        temp("{sample}.clumped.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_clumpify.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_clumpify.tsv"
    group: "bbtools"
    shell:
        """
            clumpify.sh -Xmx104g -eoom -da in1={input.raw_fwd} in2={input.raw_rev} out={output} dedupe optical 2>&1 | tee {log}
        """

rule run_bbmap_filter_by_tile:
    input:
        rules.run_bbmap_clumpify.output
    output:
        temp("{sample}.filtered_by_tile.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_filter_by_tile.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_filter_by_tile.tsv"
    group: "bbtools"
    shell:
        """
            filterbytile.sh -eoom -Xmx104g -da in={input} out={output} 2>&1 | tee {log}
        """

rule run_bbmap_bbduk_remove_adapters:
    input:
        rules.run_bbmap_filter_by_tile.output
    output:
        temp("{sample}.trimmed.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_bbduk_remove_adapters.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_bbduk_remove_adapters.tsv"
    group: "bbtools"
    shell:
        """
            bbduk.sh -Xmx104g -eoom -da in={input} out={output} ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=50 ref=adapters ftm=5 ordered 2>&1 | tee {log}
        """

rule run_bbmap_bbduk_remove_artefacts:
    input:
        rules.run_bbmap_bbduk_remove_adapters.output
    output:
        temp("{sample}.filtered.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.run_bbmap_bbduk_remove_artefacts.log"
    benchmark:
        "benchmarks/{sample}.run_bbmap_bbduk_remove_artefacts.tsv"
    group: "bbtools"
    shell:
        """
            bbduk.sh -Xmx104g -eoom -da in={input} out={output} k=31 ref=artifacts,phix ordered cardinality 2>&1 | tee {log}
        """

rule tadpole:
    input:
        rules.run_bbmap_bbduk_remove_artefacts.output
    output:
        temp("{sample}.filtered.ec.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.tadpole.log"
    benchmark:
        "benchmarks/{sample}.tadpole.tsv"
    group: "bbtools"
    shell:
        """
            tadpole.sh -Xmx104g -eoom -da in={input} out={output} mode=correct k=62 ecc ecco merge=t 2>&1 | tee {log}
        """

rule remove_human:
    input:
        rules.tadpole.output,
        rules.setup_human.output
    output:
        temp("{sample}.filtered.ec.no.hcd.fq.gz")
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.remove_human.log"
    benchmark:
        "benchmarks/{sample}.remove_human.tsv"
    group: "bbtools"
    shell:
        """
            bbmap.sh -Xmx104g -eoom -da minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
            path=../scratch qtrim=rl trimq=10 untrim -Xmx23g in={input[0]} outu={output} 2>&1 | tee {log}
        """

rule split_reads:
    input:
        rules.remove_human.output
    output:
        fwd="samples/{sample}/reads/fwd.ec.hq.fq.gz",
        rev="samples/{sample}/reads/rev.ec.hq.fq.gz"
    threads: 16
    conda:
        "../envs/conda_qc_reads.yml"
    log:
        "logs/{sample}.split_reads.log"
    benchmark:
        "benchmarks/{sample}.split_reads.tsv"
    params:
        sample_seed=42,
        subsample_count=0
    group: "bbtools"
    shell:
        """
            reformat.sh -Xmx104g -eoom -da in={input} out1={output[0]} out2={output[1]} samplereadstarget={params.subsample_count} 2>&1 | tee {log}
        """
