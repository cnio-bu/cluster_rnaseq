## Deal with optional rules
downsampling = config["parameters"]["downsampling"]["enabled"]
two_UMIs = get_params('umi_processing','umi_pattern2')

# Let stablish where are the input reads for trimming.
if UMIs:
    dir_in_trim = "umi_extract"
else:
    if downsampling: 
        dir_in_trim = "downsampled"
    else:
        dir_in_trim = "reads"


# Let stablish where are the input reads for umi extraction.
if UMIs:
    if downsampling: 
        dir_in_umi = "downsampled"
    else:
        dir_in_umi = "reads"


## Let stablish rule order whether data is single end or paired end
if downsampling:
    if single_end:
        ruleorder: downsample_single_end > downsample_paired_end
    else:
        ruleorder: downsample_paired_end > downsample_single_end

if UMIs:
    if single_end:
        ruleorder: umi_extract_single_end > umi_extract_paired_end > umi_extract_paired_end_two_umis
    elif two_UMIs != "":
        ruleorder: umi_extract_paired_end_two_umis > umi_extract_paired_end > umi_extract_single_end
    else: 
        ruleorder: umi_extract_paired_end > umi_extract_paired_end_two_umis > umi_extract_single_end

if single_end:
    ruleorder: trim_adapters_single_end > trim_adapters_paired_end
else:
    ruleorder: trim_adapters_paired_end > trim_adapters_single_end

def get_raw_fastq(wildcards,strand=1):
    if is_multi_lane(wildcards.sample):
        return f'{OUTDIR}/reads/{wildcards.sample}_R{strand}.fastq.gz'
    else:
        units.loc[wildcards.sample, "fq" + strand]


def get_R1_fragments(wildcards):
    return units.loc[wildcards.sample, 'fq1'].to_list()


def get_R2_fragments(wildcards):
    return units.loc[wildcards.sample, 'fq2'].to_list()

    
rule concat_R1_reads:
    input:get_R1_fragments
    output:
        concat_read = OUTDIR + f'/reads/{{sample}}_R1.fastq.gz'
    threads:
        get_resource('concat', 'threads')
    resources:
        mem_mb=get_resource('concat', 'mem_mb'),
        walltime=get_resource('concat', 'walltime')
    shell: 'cat {input} > {output}'


rule concat_R2_reads:
    input: get_R2_fragments
    output:
        concat_read = OUTDIR + f'/reads/{{sample}}_R2.fastq.gz'
    threads:
        get_resource('concat', 'threads')
    resources:
        mem_mb=get_resource('concat', 'mem_mb'),
        walltime=get_resource('concat', 'walltime')
    shell: 'cat {input} > {output}'


if UMIs:
    rule umi_extract_single_end:
        input:
            [OUTDIR + '/' + dir_in_umi + '/{sample}_R1.fastq.gz']
        output:
            [OUTDIR + '/umi_extract/{sample}_R1.fastq.gz']
        conda:
            "../envs/umi-tools.yaml"
        threads:
            get_resource('umi_extract_single_end', 'threads')
        resources:
            mem_mb=get_resource('umi_extract_single_end', 'mem_mb'),
            walltime=get_resource('umi_extract_single_end', 'walltime')
        params:
            get_params('umi_processing','umi_pattern1')
        log:
            f"{LOGDIR}/umi_extract/{{sample}}.log"
        shell:"""
        umi_tools extract --stdin={input} --bc-pattern={params} --log={log} --stdout={output}
    """   

    rule umi_extract_paired_end:
        input:
            f1=OUTDIR + '/' + dir_in_umi + '/{sample}_R1.fastq.gz',
            f2=OUTDIR + '/' + dir_in_umi + '/{sample}_R2.fastq.gz'
        output:
            umi_extract1=OUTDIR + '/umi_extract/{sample}_R1.fastq.gz',
            umi_extract2=OUTDIR + '/umi_extract/{sample}_R2.fastq.gz'
        conda:
            "../envs/umi-tools.yaml"
        threads:
            get_resource('umi_extract_paired_end', 'threads')
        resources:
            mem_mb=get_resource('umi_extract_paired_end', 'mem_mb'),
            walltime=get_resource('umi_extract_paired_end', 'walltime')
        params:
            get_params('umi_processing','umi_pattern1')
        log:
            f"{LOGDIR}/umi_extract/{{sample}}.log"
        shell:"""
        umi_tools extract -I {input.f1} --bc-pattern={params} --read2-in={input.f2} --log={log} \
        --stdout={output.umi_extract1} --read2-out={output.umi_extract2}
    """

    rule umi_extract_paired_end_two_umis:
        input:
            f1=OUTDIR + '/' + dir_in_umi + '/{sample}_R1.fastq.gz',
            f2=OUTDIR + '/' + dir_in_umi + '/{sample}_R2.fastq.gz'
        output:
            umi_extract1=OUTDIR + '/umi_extract/{sample}_R1.fastq.gz',
            umi_extract2=OUTDIR + '/umi_extract/{sample}_R2.fastq.gz'
        conda:
            "../envs/umi-tools.yaml"
        threads:
            get_resource('umi_extract_paired_end', 'threads')
        resources:
            mem_mb=get_resource('umi_extract_paired_end', 'mem_mb'),
            walltime=get_resource('umi_extract_paired_end', 'walltime')
        params:
            pattern1=get_params('umi_processing','umi_pattern1'),
            pattern2=get_params('umi_processing','umi_pattern2')
        log:
            f"{LOGDIR}/umi_extract/{{sample}}.log"
        shell:"""
        umi_tools extract -I {input.f1} --bc-pattern={params.pattern1} --bc-pattern2={params.pattern2} \
        --read2-in={input.f2} --log={log} --stdout={output.umi_extract1} --read2-out={output.umi_extract2}
    """


       
rule trim_adapters_single_end:
    input:
        sample=[OUTDIR + '/' + dir_in_trim + '/{sample}_R1.fastq.gz']
    output:
        trimmed=OUTDIR + '/trimmed/{sample}/{sample}_R1.fastq.gz',
        singleton=OUTDIR + '/trimmed/{sample}/{sample}.single.fastq.gz',
        discarded=OUTDIR + '/trimmed/{sample}/{sample}.discarded.fastq.gz',
        stats=OUTDIR + '/trimmed/{sample}/{sample}.stats.txt'
    threads:
        get_resource('trim_adapters_single_end', 'threads')
    resources:
        mem_mb=get_resource('trim_adapters_single_end', 'mem_mb'),
        walltime=get_resource('trim_adapters_single_end', 'walltime')
    params:
        adapters='ref=' + get_params('trimming','adapters'),
        extra=get_params('trimming', 'extra')
    log:
        f"{LOGDIR}/trim_adapters_single_end/{{sample}}.log",
    wrapper:
        "v1.23.5/bio/bbtools/bbduk"


rule trim_adapters_paired_end:
    input:
        sample=expand(OUTDIR + '/' + dir_in_trim + '/{{sample}}_R{strand}.fastq.gz', strand=[1,2])
    output:
        trimmed=expand(OUTDIR + '/trimmed/{{sample}}/{{sample}}_R{strand}.fastq.gz', strand=[1,2]),
        singleton=OUTDIR + '/trimmed/{sample}/{sample}.single.fastq.gz',
        discarded=OUTDIR + '/trimmed/{sample}/{sample}.discarded.fastq.gz',
        stats=OUTDIR + '/trimmed/{sample}/{sample}.stats.txt'
    threads:
        get_resource('trim_adapters_paired_end', 'threads')
    resources:
        mem_mb=get_resource('trim_adapters_paired_end', 'mem_mb'),
        walltime=get_resource('trim_adapters_paired_end', 'walltime')
    params:
        adapters='ref=' + get_params('trimming','adapters'),
        extra=get_params('trimming', 'extra') + ' tpe tbo'
    log:
        f"{LOGDIR}/trim_adapters_paired_end/{{sample}}.log",
    wrapper:
        "v1.23.5/bio/bbtools/bbduk"


if downsampling:
    rule downsample_single_end:
        input:
            lambda wildcards: get_raw_fastq(wildcards, strand=1)
        output:
            OUTDIR + '/downsampled/{sample}_R1.fastq.gz'
        threads:
            get_resource('downsample_single_end', 'threads')
        resources:
            mem_mb=get_resource('downsample_single_end', 'mem_mb'),
            walltime=get_resource('downsample_single_end', 'walltime')
        params:
            n=get_params('downsampling', 'n'),
            seed=get_params('downsampling', 'seed')
        log:
            f"{LOGDIR}/downsampled/{{sample}}.log"
        wrapper:
            "0.74.0/bio/seqtk/subsample/se"


if downsampling:
    rule downsample_paired_end:
        input:
            f1=lambda wildcards: get_raw_fastq(wildcards, strand=1),
            f2=lambda wildcards: get_raw_fastq(wildcards, strand=2)
        output:
            f1=OUTDIR + '/downsampled/{sample}_R1.fastq.gz',
            f2=OUTDIR + '/downsampled/{sample}_R2.fastq.gz'
        threads:
            get_resource('downsample_paired_end', 'threads')
        resources:
            mem_mb=get_resource('downsample_paired_end', 'mem_mb'),
            walltime=get_resource('downsample_paired_end', 'walltime')
        params:
            n=get_params('downsampling', 'n'),
            seed=get_params('downsampling', 'seed')
        log:
            f"{LOGDIR}/downsampled/{{sample}}.log"
        wrapper:
            "0.74.0/bio/seqtk/subsample/pe"
