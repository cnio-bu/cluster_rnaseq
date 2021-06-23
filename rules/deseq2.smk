rule htseq_count_matrix:
    input:
        quant=expand(f"{OUTDIR}/quant/{chosen_aligner}/{{sample}}.tab", \
                     sample=samples['sample'])
    output:
        counts=f"{OUTDIR}/deseq2/{chosen_aligner}/counts.tsv"
    threads:
        get_resource("htseq_count_matrix", "threads")
    resources:
        mem=get_resource("htseq_count_matrix", "mem"),
        walltime=get_resource("htseq_count_matrix", "walltime")
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/htseq_count_matrix.log"
    script:
        "../scripts/htseq_count_matrix.R"

rule deseq2_init:
    input:
        counts=f"{OUTDIR}/deseq2/{chosen_aligner}/counts.tsv"
    output:
        dds=f"{OUTDIR}/deseq2/{chosen_aligner}/dds.rds",
        normalized_counts=f"{OUTDIR}/deseq2/{chosen_aligner}/normalizedcounts.tsv",
        vst=f"{OUTDIR}/deseq2/{chosen_aligner}/dds_vst.rds"
    threads:
        get_resource("deseq2_init", "threads")
    resources:
        mem=get_resource("deseq2_init", "mem"),
        walltime=get_resource("deseq2_init", "walltime")
    params:
        samples=config['samples'],
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        ref_levels=ref_levels,
        design=config['parameters']['deseq2']['design']
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/deseq2_init.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_init.R"

rule deseq2_diffexp:
    input:
        dds=f"{OUTDIR}/deseq2/{chosen_aligner}/dds.rds"
    output:
        xlsx=f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_diffexp.xlsx",
        tsv=f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_diffexp.tsv"
    threads:
        get_resource("deseq2_diffexp", "threads")
    resources:
        mem=get_resource("deseq2_diffexp", "mem"),
        walltime=get_resource("deseq2_diffexp", "walltime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{{contrast}}/deseq2_diffexp.log"
    conda:
        "../envs/deseq2.yaml"
    script: 
        "../scripts/deseq2_diffexp.R"