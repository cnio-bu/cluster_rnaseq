rule pca:
    input:
        vst=f"{OUTDIR}/deseq2/{deseq_path}/dds_vst{{fsuffix}}.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}." + "pdf",
        png=f"{OUTDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/{{ALLcontrast}}_pca{{fsuffix,.*}}." + "png"
    threads:
        get_resource("pca", "threads")
    resources:
        mem_mb=get_resource("pca", "mem_mb"),
        runtime=get_resource("pca", "runtime")
    params:
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        design=config['parameters']['deseq2']['design'],
        ref_levels=ref_levels.to_list()
    log: f"{LOGDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/pca{{fsuffix}}.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/pca.R"

rule ma:
    input:
        dds=rules.deseq2_init.output.dds
    output:
        pdf=f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_MAplot.pdf",
        png=f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_MAplot.png"
    threads:
        get_resource("ma", "threads")
    resources:
        mem_mb=get_resource("ma", "mem_mb"),
        runtime=get_resource("ma", "runtime")
    params:
        condition=var_interest,
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: f"{LOGDIR}/deseq2/{deseq_path}/{{contrast}}/plots/MAplot.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/MAplot.R"

rule distance:
    input:
        vst=f"{OUTDIR}/deseq2/{deseq_path}/dds_vst{{fsuffix}}.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}." + "pdf",
        png=f"{OUTDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/{{ALLcontrast}}_dist{{fsuffix,.*}}." + "png"
    threads:
        get_resource("distance", "threads")
    resources:
        mem_mb=get_resource("distance", "mem_mb"),
        runtime=get_resource("distance", "runtime")
    params:
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        ref_levels=ref_levels.to_list()
    log: f"{LOGDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/dist{{fsuffix}}.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/distances.R"

rule expression_heatmap:
    input:
        normalized_counts=rules.deseq2_init.output.normalized_counts,
        diffexp=rules.deseq2_diffexp.output.tsv
    output:
        pdf=f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_topbottomDEgenes.pdf",
        png=f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_topbottomDEgenes.png"
    threads:
        get_resource("expression_heatmap", "threads")
    resources:
        mem_mb=get_resource("expression_heatmap", "mem_mb"),
        runtime=get_resource("expression_heatmap", "runtime")
    params:
        levels=lambda wildcards: contrasts[wildcards.contrast],
        designmatrix=config['parameters']['deseq2']['designmatrix'],
        ref_levels=ref_levels.to_list()
    log: f"{LOGDIR}/deseq2/{deseq_path}/{{contrast}}/plots/topbottomDEgenes.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/topbottomDEgenes.R"