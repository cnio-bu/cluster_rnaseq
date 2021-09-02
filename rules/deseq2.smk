def get_init_input():

    if chosen_aligner == "salmon":
        init_input = f"{OUTDIR}/deseq2/{chosen_aligner}/counts.tsv"
    else:
        init_input = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/counts.tsv"
    return init_input


def get_init_output():

    if chosen_aligner == "salmon":
        dds = f"{OUTDIR}/deseq2/{chosen_aligner}/dds.rds"
        normalized_counts = f"{OUTDIR}/deseq2/{chosen_aligner}/normalizedcounts.tsv"
        vst = expand(f"{OUTDIR}/deseq2/{chosen_aligner}/dds_vst{{fsuffix}}.rds", fsuffix = filesuffix)
    else:
        dds = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds.rds"
        normalized_counts = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/normalizedcounts.tsv"
        vst = expand(f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds_vst{{fsuffix}}.rds", fsuffix = filesuffix)
    return dds, normalized_counts, vst


def get_init_log():

    if chosen_aligner == "salmon":
        init_log = f"{LOGDIR}/deseq2/{chosen_aligner}/deseq2_init.log"
    else:
        init_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/deseq2_init.log"
    return init_log


def get_diffexp_output():

    if chosen_aligner == "salmon":
        xlsx = f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_diffexp.xlsx",
        tsv = f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_diffexp.tsv"
    else:
        xlsx = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/{{contrast}}_diffexp.xlsx",
        tsv = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/{{contrast}}_diffexp.tsv"
    return xlsx, tsv


def get_diffexp_log():

    if chosen_aligner == "salmon":
        diffexp_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{{contrast}}/deseq2_diffexp.log"
    else:
        diffexp_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/deseq2_diffexp.log"
    return diffexp_log


rule deseq2_init:
    input:
        counts = get_init_input()
    output:
        dds = get_init_output()[0],
        normalized_counts = get_init_output()[1],
        vst = get_init_output()[2]
    threads:
        get_resource("deseq2_init", "threads")
    resources:
        mem=get_resource("deseq2_init", "mem"),
        walltime=get_resource("deseq2_init", "walltime")
    params:
        samples=config['samples'],
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        ref_levels=ref_levels.to_list(),
        design=config['parameters']['deseq2']['design'],
        batch=batch
    log: get_init_log()
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_init.R"


rule deseq2_diffexp:
    input:
        dds = rules.deseq2_init.output.dds
    output:
        xlsx = get_diffexp_output()[0],
        tsv = get_diffexp_output()[1]
    threads:
        get_resource("deseq2_diffexp", "threads")
    resources:
        mem=get_resource("deseq2_diffexp", "mem"),
        walltime=get_resource("deseq2_diffexp", "walltime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: get_diffexp_log()
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_diffexp.R"