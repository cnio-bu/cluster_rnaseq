rule deseq2_init:
    input:
        counts = f"{OUTDIR}/deseq2/{deseq_path}/counts.tsv"
    output:
        dds = f"{OUTDIR}/deseq2/{deseq_path}/dds.rds",
        normalized_counts = f"{OUTDIR}/deseq2/{deseq_path}/normalizedcounts.tsv",
        vst = expand(f"{OUTDIR}/deseq2/{deseq_path}/dds_vst{{fsuffix}}.rds", fsuffix = filesuffix)
    threads:
        get_resource("deseq2_init", "threads")
    resources:
        mem_mb=get_resource("deseq2_init", "mem_mb"),
        runtime=get_resource("deseq2_init", "runtime")
    params:
        samples=config['samples'],
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        ref_levels=ref_levels.to_list(),
        design=config['parameters']['deseq2']['design'],
        batch=batch
    log: f"{LOGDIR}/deseq2/{deseq_path}/deseq2_init.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_init.R"


rule deseq2_diffexp:
    input:
        dds = rules.deseq2_init.output.dds
    output:
        xlsx = f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.xlsx",
        tsv = f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.tsv",
        xlsx_lfcShrink = f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_lfcShrink_diffexp.xlsx",
        tsv_lfcShrink = f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_lfcShrink_diffexp.tsv"
    threads:
        get_resource("deseq2_diffexp", "threads")
    resources:
        mem_mb=get_resource("deseq2_diffexp", "mem_mb"),
        runtime=get_resource("deseq2_diffexp", "runtime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: f"{LOGDIR}/deseq2/{deseq_path}/{{contrast}}/deseq2_diffexp.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_diffexp.R"