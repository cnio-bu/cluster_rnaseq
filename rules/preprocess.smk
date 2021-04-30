## Let snakemake use paired reads whenever possible
ruleorder: trim_adapters_paired_end > trim_adapters_single_end


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), 'fq2']).all()


def is_multi_lane(sample, unit):
    return not pd.isnull(units.loc[(sample,unit), 'lane']).all()


def get_single_raw_fastq(wildcards):

    if is_multi_lane(wildcards.sample, wildcards.unit):
        return f'{OUTDIR}/reads/{wildcards.sample}_{wildcards.unit}_R1_concat.fastq.gz'
    else:
        return units.loc[(wildcards.sample, wildcards.unit), "fq1"]


def get_paired_raw_fastq(wildcards):
    
    if is_multi_lane(wildcards.sample, wildcards.unit):
        return expand(OUTDIR + '/reads/{sample}_{unit}_R{strand}_concat.fastq.gz', strand=[1,2], **wildcards)
    else:
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


def get_R1_fragments(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), 'fq1'].to_list()


def get_R2_fragments(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), 'fq2'].to_list()


rule concat_R1_reads:
	input:get_R1_fragments
	output:
		concat_read = OUTDIR + f'/reads/{{sample}}_{{unit}}_R1_concat.fastq.gz'
	shell: 'cat {input} > {output}'


rule concat_R2_reads:
	input: get_R2_fragments
	output:
		concat_read = OUTDIR + f'/reads/{{sample}}_{{unit}}_R2_concat.fastq.gz'
	shell: 'cat {input} > {output}'


rule trim_adapters_single_end:
    input:
        get_single_raw_fastq,
    output:
        fastq=OUTDIR + "/trimmed/{sample}_{unit}_R1.fastq.gz",
        qc=OUTDIR + "/qc/trimming/{sample}_{unit}.qc.txt"
    threads:
        get_resource('trim_adapters_single_end', 'threads')
    resources:
        mem=get_resource('trim_adapters_single_end', 'mem'),
        walltime=get_resource('trim_adapters_single_end', 'walltime')
    params:
        adapters=f"-a {{get_params(trimming,'adapters')}}",
        extra=get_params('trimming', 'extra')
    log:
        "logs/trim_adapters_single_end/{sample}-{unit}.log",
    wrapper:
        '0.74.0/bio/cutadapt/se'


rule trim_adapters_paired_end:
    input:
        get_paired_raw_fastq,
    output:
        fastq1=OUTDIR + "/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2=OUTDIR + "/trimmed/{sample}_{unit}_R2.fastq.gz",
        qc=OUTDIR + "/qc/trimming/{sample}_{unit}.qc.txt"
    threads:
        get_resource('trim_adapters_single_end', 'threads')
    resources:
        mem=get_resource('trim_adapters_single_end', 'mem'),
        walltime=get_resource('trim_adapters_single_end', 'walltime')
    params:
        adapters=f"-a {{get_params(trimming,'adapters')}}",
        extra=get_params('trimming', 'extra')
    log:
        "logs/trim_adapters_paired_end/{sample}-{unit}.log",
    wrapper:
        '0.74.0/bio/cutadapt/pe'