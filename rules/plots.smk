rule pca:
    input:
        vst=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/dds_vst{{fsuffix}}.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/{{ALLcontrast}}_pca{{fsuffix,.*}}.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/{{ALLcontrast}}_pca{{fsuffix,.*}}.png"
    #threads:
    #    get_resource("deseq2_diffexp", "threads")
    #resources:
    #    mem=get_resource("deseq2_diffexp", "mem"),
    #    walltime=get_resource("deseq2_diffexp", "walltime")
    params:
        levels=lambda wildcards: allSamples[wildcards.ALLcontrast],
        design=config['parameters']['deseq2']['design']
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{chosen_quantifier}/{{ALLcontrast}}/pca{{fsuffix,.*}}.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/pca.R"
