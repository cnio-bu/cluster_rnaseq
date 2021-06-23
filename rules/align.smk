## Let snakemake use paired reads whenever possible
ruleorder: salmon_quant_paired > salmon_quant_se
ruleorder: star_align_paired > star_align_se

def get_hisat_reads(wildcards):
    if is_single_end(wildcards.sample):
        return f"{OUTDIR}/trimmed/{wildcards.sample}_R1.fastq.gz"
    else:
        return expand(f"{OUTDIR}/trimmed/{{sample}}_R{{strand}}.fastq.gz", strand=[1,2], **wildcards)


### SALMON ###
rule salmon_quant_se:
    input:
        salmon_index=rules.salmon_index.output,
        se_reads=rules.trim_adapters_single_end.output.fastq
    output:
        quant=directory(OUTDIR + '/quant/salmon/{sample}/')
    threads:
        get_resource('salmon_quant', 'threads')
    resources:
        mem=get_resource('salmon_quant', 'mem'),
        walltime=get_resource('salmon_quant', 'walltime')
    params:
        libtype=get_params('salmon', 'libtype')
    log:
        f"{LOGDIR}/salmon_quant/{{sample}}.log"
    conda:
        '../envs/aligners.yaml'
    shell: 'salmon quant -i {input.salmon_index} -l {params.libtype} -r {input.se_reads} --validateMappings -o {output.quant} --threads {threads}'


rule salmon_quant_paired:
    input:
        salmon_index=rules.salmon_index.output,
        r1_reads=rules.trim_adapters_paired_end.output.fastq1,
        r2_reads=rules.trim_adapters_paired_end.output.fastq2
    output:
        quant=directory(OUTDIR + '/quant/salmon/{sample}/')
    threads:
        get_resource('salmon_quant', 'threads')
    resources:
        mem=get_resource('salmon_quant', 'mem'),
        walltime=get_resource('salmon_quant', 'walltime')
    params:
        libtype=get_params('salmon', 'libtype')
    log:
        f"{LOGDIR}/salmon_quant/{{sample}}.log"
    conda:
        '../envs/aligners.yaml'
    shell: 'salmon quant -i {input.salmon_index} -l {params.libtype} -1 {input.r1_reads} -2 {input.r2_reads} --validateMappings -o {output.quant} --threads {threads}'


### STAR ###
rule star_align_se:
    input:
        fq1=rules.trim_adapters_single_end.output.fastq
    output:
        aligned=OUTDIR + '/mapped/star/{sample}/Aligned.sortedByCoord.out.bam'
    threads:
        get_resource('star_align', 'threads')
    resources:
        mem=get_resource('star_align', 'mem'),
        walltime=get_resource('star_align', 'walltime')
    params:
        # path to STAR reference genome index
        index=config["ref"]["star"]["star_index"],
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate"
    log:
        f"{LOGDIR}/star/{{sample}}.log"
    wrapper:
        "0.74.0/bio/star/align"


rule star_align_paired:
    input:
        fq1=rules.trim_adapters_paired_end.output.fastq1,
        fq2=rules.trim_adapters_paired_end.output.fastq2
    output:
        aligned=OUTDIR + '/mapped/star/{sample}/Aligned.sortedByCoord.out.bam'
    threads:
        get_resource("star_align", "threads")
    resources:
        mem=get_resource('star_align', 'mem'),
        walltime=get_resource('star_align', 'walltime')
    params:
        # path to STAR reference genome index
        index=config["ref"]["star"]["star_index"],
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate"
    log:
        f"{LOGDIR}/star/{{sample}}.log"
    wrapper:
        "0.74.0/bio/star/align"


## HISAT-2 ##
rule hisat2_align:
    input:
        reads=get_hisat_reads
    output:
        aligned=OUTDIR + "/mapped/hisat2/{sample}.bam"
    log:
        f"{LOGDIR}/hisat2_align/{{sample}}.log"
    params:
        extra="--new-summary",
        idx=config["ref"]["hisat2"]["hisat2_index"] + "/hisat2_index"
    threads:
        get_resource("hisat2_align", "threads")
    resources:
        mem=get_resource('hisat2_align', 'mem'),
        walltime=get_resource('hisat2_align', 'walltime')
    wrapper:
        "0.74.0/bio/hisat2/align"


rule hisat2_sort:
    input:
        aligned=OUTDIR + "/mapped/hisat2/{sample}.bam"
    output:
        sortedCoord=OUTDIR + "/mapped/hisat2/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        f"{LOGDIR}/hisat2_sort/{{sample}}_sorted.log"
    threads:
        get_resource("hisat2_sort", "threads")
    resources:
        mem=get_resource('hisat2_sort', 'mem'),
        walltime=get_resource('hisat2_sort', 'walltime')
    conda:
        '../envs/aligners.yaml'
    shell:
        'samtools sort -@ {threads} -o {output.sortedCoord} {input.aligned}'