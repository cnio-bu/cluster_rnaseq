def get_vst_input():

    if chosen_aligner == "salmon":
        vst_input = f"{OUTDIR}/deseq2/{chosen_aligner}/dds_vst{{fsuffix}}.rds"
    else:
        vst_input = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds_vst{{fsuffix}}.rds"
    return vst_input

def get_pca_output(ext):

    if chosen_aligner == "salmon":
        pca_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}." + ext
    else:
        pca_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}." + ext
    return pca_output

def get_pca_log():

    if chosen_aligner == "salmon":
        pca_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{{ALLcontrast}}/plots/pca{{fsuffix}}.log"
    else:
        pca_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/pca{{fsuffix}}.log"
    return pca_log

def get_distance_output(ext):

    if chosen_aligner == "salmon":
        distance_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}." + ext
    else:
        distance_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}." + ext
    return distance_output

def get_distance_log():

    if chosen_aligner == "salmon":
        distance_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{{ALLcontrast}}/plots/dist{{fsuffix}}.log"
    else:
        distance_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/dist{{fsuffix}}.log"
    return distance_log

def get_contrast_output(plot, ext):

    if chosen_aligner == "salmon":
        ma_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/plots/{{contrast}}_" + plot + "." + ext
    else:
        ma_output = f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/{{contrast}}_" + plot + "." + ext
    return ma_output

def get_contrast_log(plot):

    if chosen_aligner == "salmon":
        contrast_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{{contrast}}/plots/" + plot + ".log"
    else:
        contrast_log = f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/" + plot + ".log"
    return contrast_log


rule pca:
    input:
        vst=get_vst_input()
    output:
        pdf=get_pca_output(ext="pdf"),
        png=get_pca_output(ext="png")
    threads:
        get_resource("pca", "threads")
    resources:
        mem=get_resource("pca", "mem"),
        walltime=get_resource("pca", "walltime")
    params:
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        design=config['parameters']['deseq2']['design']
    log: get_pca_log()
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/pca.R"

rule ma:
    input:
        dds=rules.deseq2_init.output.dds
    output:
        pdf=get_contrast_output(plot="MAplot", ext="pdf"),
        png=get_contrast_output(plot="MAplot", ext="png"),
    threads:
        get_resource("ma", "threads")
    resources:
        mem=get_resource("ma", "mem"),
        walltime=get_resource("ma", "walltime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: get_contrast_log(plot="MAplot")
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/MAplot.R"

rule distance:
    input:
        vst=get_vst_input()
    output:
        pdf=get_distance_output(ext="pdf"),
        png=get_distance_output(ext="png")
    threads:
        get_resource("distance", "threads")
    resources:
        mem=get_resource("distance", "mem"),
        walltime=get_resource("distance", "walltime")
    params:
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast]
    log: get_distance_log()
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/distances.R"

rule expression_heatmap:
    input:
        normalized_counts=rules.deseq2_init.output.normalized_counts,
        diffexp=rules.deseq2_diffexp.output.tsv
    output:
        pdf=get_contrast_output(plot="topbottomDEgenes", ext="pdf"),
        png=get_contrast_output(plot="topbottomDEgenes", ext="png")
    threads:
        get_resource("expression_heatmap", "threads")
    resources:
        mem=get_resource("expression_heatmap", "mem"),
        walltime=get_resource("expression_heatmap", "walltime")
    params:
        levels=lambda wildcards: contrasts[wildcards.contrast],
        designmatrix=config['parameters']['deseq2']['designmatrix']
    log: get_contrast_log(plot="topbottomDEgenes")
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/topbottomDEgenes.R"