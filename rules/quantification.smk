def featurecounts_args(sample):
    pars = ""
    single_end = is_single_end(sample)
    if single_end == False:
        pars = "-p"
    return pars


###  BAM INDEXING ###
rule bam_indexing:
    input:
        aligned=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        bai_index=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"
    log:
        f"{LOGDIR}/bam_indexing/{{sample}}.log"
    threads:
        get_resource("bam_indexing", "threads")
    resources:
        mem=get_resource('bam_indexing', 'mem'),
        walltime=get_resource('bam_indexing', 'walltime')
    conda:
        '../envs/aligners.yaml'
    shell:
        'samtools index -@ {threads} {input.aligned}'


### HTSEQ COUNT ###
rule htseq_count:
    input:
        bam_file=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam",
        bai_index=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"
    output:
        quant=f"{OUTDIR}/quant/{chosen_aligner}/htseq/{{sample}}.tab"
    threads:
        get_resource('htseq_count', 'threads')
    resources:
        mem=get_resource('htseq_count', 'mem'),
        walltime=get_resource('htseq_count', 'walltime')
    params:
        annotation= lambda x: config['ref'][chosen_aligner]['annotation'] if chosen_aligner != 'salmon' else '',
        extra='-f bam -r pos' + config['parameters']['htseq-count']['extra'],
        mode = config['parameters']['htseq-count']['mode'],
        strandedness = config['parameters']['htseq-count']['strandedness']
    log:
        f"{LOGDIR}/htseq_count/{{sample}}.log"
    conda:
        '../envs/cuantification.yaml'
    shell: 'htseq-count {params.extra} {params.mode} {params.strandedness} {input.bam_file} {params.annotation} > {output.quant} 2> {log}'


rule htseq_count_matrix:
    input:
        quant=expand(f"{OUTDIR}/quant/{chosen_aligner}/htseq/{{sample}}.tab", \
                     sample=samples['sample'])
    output:
        counts=f"{OUTDIR}/deseq2/{chosen_aligner}/htseq/counts.tsv"
    threads:
        get_resource("htseq_count_matrix", "threads")
    resources:
        mem=get_resource("htseq_count_matrix", "mem"),
        walltime=get_resource("htseq_count_matrix", "walltime")
    log: f"{LOGDIR}/deseq2/{chosen_aligner}/htseq_count_matrix.log"
    conda:
        '../envs/deseq2.yaml'
    script:
        "../scripts/htseq_count_matrix.R"


### FEATURECOUNTS ###
rule featurecounts:
    input:
        bam_file= f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        quant=f"{OUTDIR}/quant/{chosen_aligner}/featureCounts/{{sample}}.tab"
    threads:
        get_resource('featureCounts', 'threads')
    resources:
        mem=get_resource('featureCounts', 'mem'),
        walltime=get_resource('featureCounts', 'walltime')
    params:
        args= lambda wc: featurecounts_args(wc.sample),
        extra= config['parameters']['featureCounts']['extra'],
        annotation= lambda x: config['ref'][chosen_aligner]['annotation'] if chosen_aligner != 'salmon' else '',
        strandedness = config['parameters']['featureCounts']['strandedness']
    log:
        f"{LOGDIR}/featureCounts/{{sample}}.log"
    conda:
        '../envs/cuantification.yaml'
    shell: 'featureCounts {params.args} -T {threads} {params.extra} {params.strandedness} -a {params.annotation} -o {output.quant} {input.bam_file} &> {log}'


rule fcounts_count_matrix:
    input:
        expand(f"{OUTDIR}/quant/{chosen_aligner}/featureCounts/{{sample}}.tab", sample=samples['sample'])
    output:
        counts=f"{OUTDIR}/deseq2/{chosen_aligner}/featureCounts/counts.tsv"
    threads:
        get_resource('fcounts_count_matrix', 'threads')
    resources:
        mem=get_resource('fcounts_count_matrix', 'mem'),
        walltime=get_resource('fcounts_count_matrix', 'walltime')
    params:
        samples=config['samples'],
    script:
        "../scripts/fcounts_count_matrix.py"


### SALMON MATRIX ###
rule salmon_matrix_from_quants:
    input:
        quants = expand(f"{OUTDIR}/quant/salmon/{{sample}}/quant.sf",  sample=samples['sample'])
    output:
        gene_level_matrix    = f"{OUTDIR}/deseq2/salmon/counts.tsv",
        transcript_estimates = f"{OUTDIR}/deseq2/salmon/transcript_level_estimates.rds",
        metadata_cache      = temp(directory(f"{OUTDIR}/quant/salmon/metadata_cache"))

    threads:  
        get_resource('salmon_matrix_from_quants', 'threads')
    resources:
        mem=get_resource('salmon_matrix_from_quants', 'mem'),
        walltime=get_resource('salmon_matrix_from_quants', 'walltime')
    params:
        salmon_quant_directory = f"{OUTDIR}/quant/salmon",
        samples                  = config['samples']
    log:
        f"{LOGDIR}/salmon_matrix_from_quants.log"
    conda:
        '../envs/cuantification.yaml'
    script:
        '../scripts/salmon_matrix_from_quant.R'
