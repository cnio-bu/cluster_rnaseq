## Let snakemake use paired reads whenever possible
ruleorder: trim_adapters_paired_end > trim_adapters_single_end

def get_single_raw_fastq(wildcards):

    if is_multi_lane(wildcards.sample):
        return [f'{OUTDIR}/reads/{wildcards.sample}_R1_concat.fastq.gz']
    else:
        return [units.loc[wildcards.sample, "fq1"]]


def get_paired_raw_fastq(wildcards):
    
    if is_multi_lane(wildcards.sample):
        return expand(OUTDIR + '/reads/{sample}_R{strand}_concat.fastq.gz', strand=[1,2], **wildcards)
    else:
        return units.loc[wildcards.sample, ["fq1", "fq2"]].dropna().to_list()


def get_R1_fragments(wildcards):
    return units.loc[wildcards.sample, 'fq1'].to_list()


def get_R2_fragments(wildcards):
    return units.loc[wildcards.sample, 'fq2'].to_list()


rule concat_R1_reads:
	input:get_R1_fragments
	output:
		concat_read = OUTDIR + f'/reads/{{sample}}_R1_concat.fastq.gz'
	shell: 'cat {input} > {output}'


rule concat_R2_reads:
	input: get_R2_fragments
	output:
		concat_read = OUTDIR + f'/reads/{{sample}}_R2_concat.fastq.gz'
	shell: 'cat {input} > {output}'


rule trim_adapters_single_end:
    input:
        sample=get_single_raw_fastq
    output:
        trimmed=OUTDIR + '/trimmed/{sample}/{sample}_R1.fastq.gz',
        singleton=OUTDIR + '/trimmed/{sample}/{sample}.single.fastq.gz',
        discarded=OUTDIR + '/trimmed/{sample}/{sample}.discarded.fastq.gz',
        stats=OUTDIR + '/trimmed/{sample}/{sample}.stats.txt'
    threads:
        get_resource('trim_adapters_single_end', 'threads')
    resources:
        mem=get_resource('trim_adapters_single_end', 'mem'),
        walltime=get_resource('trim_adapters_single_end', 'walltime')
    params:
        adapters='ref=' + get_params('trimming','adapters'),
        extra=get_params('trimming', 'extra')
    log:
        f"{LOGDIR}/trim_adapters_single_end/{{sample}}.log",
    wrapper:
        '0.74.0/bio/bbtools/bbduk'


rule trim_adapters_paired_end:
    input:
        sample=get_paired_raw_fastq
    output:
        trimmed=expand(OUTDIR + '/trimmed/{{sample}}/{{sample}}_R{strand}.fastq.gz', strand=[1,2]),
        singleton=OUTDIR + '/trimmed/{sample}/{sample}.single.fastq.gz',
        discarded=OUTDIR + '/trimmed/{sample}/{sample}.discarded.fastq.gz',
        stats=OUTDIR + '/trimmed/{sample}/{sample}.stats.txt'
    threads:
        get_resource('trim_adapters_paired_end', 'threads')
    resources:
        mem=get_resource('trim_adapters_paired_end', 'mem'),
        walltime=get_resource('trim_adapters_paired_end', 'walltime')
    params:
        adapters='ref=' + get_params('trimming','adapters'),
        extra=get_params('trimming', 'extra') + ' tpe tbo'
    log:
        f"{LOGDIR}/trim_adapters_paired_end/{{sample}}.log",
    wrapper:
        '0.74.0/bio/bbtools/bbduk'