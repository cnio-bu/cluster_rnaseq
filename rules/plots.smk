rule pca:
    input:
        vst=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds_vst{{fsuffix}}.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}.png"
    threads:
        get_resource("pca", "threads")
    resources:
        mem=get_resource("pca", "mem"),
        walltime=get_resource("pca", "walltime")
    params:
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        design=config['parameters']['deseq2']['design']
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/pca{{fsuffix,.*}}.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/pca.R"

rule ma:
    input:
        dds=rules.deseq2_init.output.dds
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/{{contrast}}_MAplot.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/{{contrast}}_MAplot.png"
    threads:
        get_resource("ma", "threads")
    resources:
        mem=get_resource("ma", "mem"),
        walltime=get_resource("ma", "walltime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/MAplot.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/MAplot.R"

rule distance:
    input:
        vst=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds_vst{{fsuffix}}.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}.png"
    threads:
        get_resource("distance", "threads")
    resources:
        mem=get_resource("distance", "mem"),
        walltime=get_resource("distance", "walltime")
    params:
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        designmatrix=config['parameters']['deseq2']['designmatrix']
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/plots/dist{{fsuffix,.*}}.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/distances.R"

rule expression_heatmap:
    input:
        normalized_counts=rules.deseq2_init.output.normalized_counts,
        diffexp=rules.deseq2_diffexp.output.tsv
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/{{contrast}}_topbottomDEgenes.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/{{contrast}}_topbottomDEgenes.png"
    threads:
        get_resource("expression_heatmap", "threads")
    resources:
        mem=get_resource("expression_heatmap", "mem"),
        walltime=get_resource("expression_heatmap", "walltime")
    params:
        levels=lambda wildcards: contrasts[wildcards.contrast],
        designmatrix=config['parameters']['deseq2']['designmatrix']
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{contrast}}/plots/topbottomDEgenes.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/topbottomDEgenes.R"