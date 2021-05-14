
def create_list(df_row, read):
    
    return f"{OUTDIR}/qc/fastqc/{df_row.loc['sample']}.{df_row.loc['unit']}.{df_row.loc['lane']}.fq{read}_fastqc.html"


def fastqc_output(units):
    
    paths_fq1 = units.apply(create_list, read="1", axis = 1).tolist()
    paths_fq2 = units[units['fq2'].notnull()].apply(create_list, read="2", axis = 1).tolist()
    all_paths = paths_fq1+paths_fq2

    return all_paths


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
        fastq= lambda wc: units.loc(axis=0)[(wc.sample,wc.unit,wc.lane)]['fq' + wc.read],
        flag=f"{LOGDIR}/md5.checked"
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}.{{unit}}.{{lane}}.fq{{read}}_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}.{{unit}}.{{lane}}.fq{{read}}_fastqc.zip"
    threads: 
        get_resource("fastqc","threads")
    resources:
        mem=get_resource("fastqc","mem"),
        walltime=get_resource("fastqc","walltime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        f"{LOGDIR}/fastqc/{{sample}}.{{unit}}.{{lane}}.fq{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastqc/{{sample}}.{{unit}}.{{lane}}.fq{{read}}.bmk"
    wrapper:
        "0.74.0/bio/fastqc"