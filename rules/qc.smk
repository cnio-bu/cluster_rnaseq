def create_list(df_row, read):
    return f"{OUTDIR}/qc/fastqc/{df_row.loc['sample']}_{df_row.loc['lane']}_fq{read}_fastqc.zip"
   

def create_list_concat(df_row, read):
    return f"{OUTDIR}/qc_concat/fastqc/{df_row.loc['sample']}_R{read}_fastqc.zip"


def fastqc_output(units):
    paths_fq1 = units.apply(create_list, read="1", axis = 1).values.tolist()
    paths_fq2 = units[units['fq2'].notnull()].apply(create_list, read="2", axis = 1).values.tolist()
    all_paths = paths_fq1+paths_fq2
    return all_paths


def multiqc_concat_input(units):
    paths_fqc_fq1 = np.unique(units.apply(create_list_concat, read="1", axis = 1).values).tolist()
    paths_fqc_fq2 = np.unique(units[units['fq2'].notnull()].apply(create_list_concat, read="2", axis = 1).values).tolist()
    all_input = paths_fqc_fq1+paths_fqc_fq2

    if config["aligner"] == 0:
        all_input += expand(f"{OUTDIR}/mapped/star/{{samples.sample}}/Aligned.sortedByCoord.out.bam", samples=samples.itertuples())
    elif config["aligner"] == 1:
        all_input += expand(f"{OUTDIR}/quant/salmon/{{samples.sample}}/quant.sf", samples=samples.itertuples())
    elif config["aligner"] == 2:
        all_input += expand(f"{LOGDIR}/hisat2_align/{{samples.sample}}.log", samples=samples.itertuples())
    else:
        "No aligner has been specified for multiQC"
    return all_input


rule md5:
    input:
        units = config["units"]
    output:
        temp(touch(f"{LOGDIR}/md5.checked"))
    params:
        logdir=LOGDIR
    conda:
        '../envs/md5.yaml'
    script:
        "../scripts/md5.py"


rule fastqc:
    input:
        fastq= lambda wc: units.loc(axis=0)[(wc.sample,wc.lane)]['fq' + wc.read],
        flag=f"{LOGDIR}/md5.checked"
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}_{{lane}}_fq{{read}}_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}_{{lane}}_fq{{read}}_fastqc.zip"
    threads: 
        get_resource("fastqc","threads")
    resources:
        mem=get_resource("fastqc","mem"),
        walltime=get_resource("fastqc","walltime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        f"{LOGDIR}/fastqc/{{sample}}_{{lane}}_fq{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastqc/{{sample}}_{{lane}}_fq{{read}}.bmk"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        fastqc_output(units)
    output:
        report(f"{OUTDIR}/qc/multiqc_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"]
    benchmark:
        f"{LOGDIR}/multiqc.bmk"
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem=get_resource("multiqc","mem"),
        walltime=get_resource("multiqc","walltime")
    wrapper:
        "0.74.0/bio/multiqc"


# QC for concat files

rule fastqc_concat:
    input:
        f"{OUTDIR}/trimmed/{{sample}}/{{sample}}_R{{read}}.fastq.gz"
    output:
        html=f"{OUTDIR}/qc_concat/fastqc/{{sample}}_R{{read}}_fastqc.html",
        zip=f"{OUTDIR}/qc_concat/fastqc/{{sample}}_R{{read}}_fastqc.zip"
    threads: 
        get_resource("fastqc","threads")
    resources:
        mem=get_resource("fastqc","mem"),
        walltime=get_resource("fastqc","walltime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        f"{LOGDIR}/fastqc/{{sample}}_R{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastqc/{{sample}}_R{{read}}.bmk"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc_concat:
    input:
        multiqc_concat_input(units)
    output:
        report(f"{OUTDIR}/qc_concat/multiqc_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"]
    benchmark:
        f"{LOGDIR}/multiqc.bmk"
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem=get_resource("multiqc","mem"),
        walltime=get_resource("multiqc","walltime")
    wrapper:
        "0.74.0/bio/multiqc"
