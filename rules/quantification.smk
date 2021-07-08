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


rule htseq_count:
    input:
        bam_file=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam",
        bai_index=f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"
    output:
        quant=f"{OUTDIR}/quant/{chosen_aligner}/{{sample}}.tab"
    threads:
        get_resource('htseq_count', 'threads')
    resources:
        mem=get_resource('htseq_count', 'mem'),
        walltime=get_resource('htseq_count', 'walltime')
    params:
        annotation= lambda x: config['ref'][chosen_aligner]['annotation'] if chosen_aligner != 'salmon' else '',
        extra='-f bam -r pos',
        mode = config['parameters']['htseq-count']['mode'],
        strandedness = config['parameters']['htseq-count']['strandedness']
    log:
        f"{LOGDIR}/htseq_count/{{sample}}.log"
    conda:
        '../envs/cuantification.yaml'
    shell: 'htseq-count {params.extra} {params.mode} {params.strandedness} {input.bam_file} {params.annotation} > {output.quant} 2> {log}'


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
