rule pca:
    input:
        vst=f"{OUTDIR}/deseq2/{chosen_aligner}/dds_vst.rds"
    output:
        pdf=f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_pca.pdf",
        png=f"{OUTDIR}/deseq2/{chosen_aligner}/{{contrast}}/{{contrast}}_pca.png"
    #threads:
    #    get_resource("deseq2_diffexp", "threads")
    #resources:
    #    mem=get_resource("deseq2_diffexp", "mem"),
    #    walltime=get_resource("deseq2_diffexp", "walltime")
    params:
        levels=lambda wildcards: contrasts[wildcards.contrast]
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/{{contrast}}/pca.log"
    conda:
        "../envs/plots.yaml"
    script: 
        "../scripts/pca.R"