chosen_aligner = get_aligner(int(config['aligner']))

### BAM INDEXING ###
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
        annotation= config['ref'][chosen_aligner]['annotation'],
        extra='-f bam -r pos',
        mode = config['parameters']['htseq-count']['mode'],
        strandedness = config['parameters']['htseq-count']['strandedness']
    log:
        f"{LOGDIR}/htseq_count/{{sample}}.log"
    conda:
        '../envs/cuantification.yaml'
    shell: 'htseq-count {params.extra} {params.mode} {params.strandedness} {input.bam_file} {params.annotation} > {output.quant}'