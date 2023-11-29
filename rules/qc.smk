## Deal with optional rules (fastq_screen)
fastq_screen = config["parameters"]["fastq_screen"]["enabled"]

def create_list_files_fastqc(df_row, read):
    return f"{OUTDIR}/qc/fastqc_files/{df_row.loc['sample']}_{df_row.loc['lane']}_fq{read}_fastqc.zip"

def create_list_files_fastq_screen(df_row, read):
    return  f"{OUTDIR}/fastq_screen/fastq_screen_files/{df_row.loc['sample']}_{df_row.loc['lane']}_fq{read}_fastq_screen.txt"

def create_list_concat_fastqc(df_row, read):
    return f"{OUTDIR}/qc/fastqc_concat/{df_row.loc['sample']}_R{read}_fastqc.zip"

def create_list_concat_fastq_screen(df_row, read):
    return f"{OUTDIR}/fastq_screen/fastq_screen_concat/{df_row.loc['sample']}_R{read}_fastq_screen.txt"

   
def multiqc_files_input(units):
    paths_fq1_fastqc = units.apply(create_list_files_fastqc, read="1", axis = 1).values.tolist()
    paths_fq2_fastqc = units[units['fq2'].notnull()].apply(create_list_files_fastqc, read="2", axis = 1).values.tolist()
    all_paths = paths_fq1_fastqc+paths_fq2_fastqc
    # Only if Fastq_screen is performed (enabled = True), fastq_screen files will be 
    # included in the multiqc_files_report.
    if fastq_screen:
        paths_fq1_fastq_screen = units.apply(create_list_files_fastq_screen, read="1", axis = 1).values.tolist()
        paths_fq2_fastq_screen = units[units['fq2'].notnull()].apply(create_list_files_fastq_screen, read="2", axis = 1).values.tolist()
        all_paths += paths_fq1_fastq_screen+paths_fq2_fastq_screen
    return all_paths


def multiqc_concat_input(units, step):

    # steps = ["all", "index", "files_qc", "trimming", "alignment", "quantification", \
    #          "diffexp", "plots", "multiqc_all"]

    if step == "index" or step == "files_qc" or step == "multiqc_all":
        return []

    # Quality Control for concat files --> FastQC --> Always included
    paths_fqc_fq1 = np.unique(units.apply(create_list_concat_fastqc, read="1", axis = 1).values).tolist()
    paths_fqc_fq2 = np.unique(units[units['fq2'].notnull()].apply(create_list_concat_fastqc, read="2", axis = 1).values).tolist()
    mqc_input = paths_fqc_fq1+paths_fqc_fq2

    # Fasq_screen for concat files --> Fastq_screen --> Only included if the enabled parameter of fastq_screen is True.
    if fastq_screen: 
        paths_fq1_fastq_screen = np.unique(units.apply(create_list_concat_fastq_screen, read="1", axis = 1).values).tolist()
        paths_fq2_fastq_screen = np.unique(units[units['fq2'].notnull()].apply(create_list_concat_fastq_screen, read="2", axis = 1).values).tolist()
        mqc_input += paths_fq1_fastq_screen+paths_fq2_fastq_screen 

    # Trimming --> BBDUK
    mqc_input += expand(f"{OUTDIR}/trimmed/{{samples.sample}}/{{samples.sample}}.stats.txt", samples=samples.itertuples())

    if step == "trimming":
        return mqc_input
    
    # Aligment
    if chosen_aligner == "star":
        mqc_input += expand(f"{OUTDIR}/mapped/star/{{samples.sample}}/Aligned.sortedByCoord.out.bam", samples=samples.itertuples())
    elif chosen_aligner == "salmon":
        mqc_input += expand(f"{OUTDIR}/quant/salmon/{{samples.sample}}/quant.sf", samples=samples.itertuples())
    elif chosen_aligner == "hisat2":
        mqc_input += expand(f"{LOGDIR}/hisat2_align/{{samples.sample}}.log", samples=samples.itertuples())
    else:
        "No aligner has been specified for multiQC"

    if step == "alignment":
        return mqc_input
    
    # Quantification
    if chosen_aligner == "star" or chosen_aligner == "hisat2":
        if chosen_quantifier == "htseq":
            mqc_input += expand(f"{OUTDIR}/quant/{chosen_aligner}/htseq/{{samples.sample}}.tab", samples=samples.itertuples())
        elif chosen_quantifier == "featureCounts":
            mqc_input += expand(f"{OUTDIR}/quant/{chosen_aligner}/featureCounts/{{samples.sample}}.tab.summary", samples=samples.itertuples())
        else:
            "No quantifier has been specified for multiQC"

    if step == "quantification" or step == "diffexp" or step == "plots":
        return mqc_input


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


# QC for single files

rule fastqc_files:
    input:
        fastq= lambda wc: units.loc(axis=0)[(wc.sample,wc.lane)]['fq' + wc.read],
        flag=f"{LOGDIR}/md5.checked"
    output:
        html=f"{OUTDIR}/qc/fastqc_files/{{sample}}_{{lane}}_fq{{read}}_fastqc.html",
        zip_qc=f"{OUTDIR}/qc/fastqc_files/{{sample}}_{{lane}}_fq{{read}}_fastqc.zip"
    threads: 
        get_resource("fastqc","threads")
    resources:
        mem_mb=get_resource("fastqc","mem_mb"),
        runtime=get_resource("fastqc","runtime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        f"{LOGDIR}/fastqc_files/{{sample}}_{{lane}}_fq{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastqc_files/{{sample}}_{{lane}}_fq{{read}}.bmk"
    wrapper:
        "0.74.0/bio/fastqc"


# The rule Fastq_screen_indexes downloads pre-built Bowtie2 indices of commonly used genomes.
# The genome indices will be downloaded to a folder named "FastQ_Screen_Genomes",  
# where the fast_screen.conf is already placed. This config file is ready to use and 
# lists the correct paths to the newly downloaded reference genomes.

rule fastq_screen_indexes:
    output:
        "res/FastQ_Screen_Genomes/fastq_screen.conf"
    conda:
        "../envs/fastq_screen.yaml"
    threads: 
        get_resource("fastq_screen_indexes","threads")
    resources:
        mem_mb=get_resource("fastq_screen_indexes","mem_mb"),
        runtime=get_resource("fastq_screen_indexes","runtime")
    params:
        outdir=config["parameters"]["fastq_screen_indexes"]["outdir"]
    log:
        f"{LOGDIR}/fastq_screen_indexes.log"
    benchmark:
        f"{LOGDIR}/fastq_screen_indexes.bmk"
    shell:"""
        fastq_screen --threads {threads} --get_genomes --outdir {params.outdir}/ &> {log}
    """



# Fastq_screen for single files.

rule fastq_screen_files:
    input: 
        fastq= lambda wc: units.loc(axis=0)[(wc.sample,wc.lane)]['fq' + wc.read],
        conf="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["parameters"]["fastq_screen_indexes"]["outdir"])
    output:
        txt=f"{OUTDIR}/fastq_screen/fastq_screen_files/{{sample}}_{{lane}}_fq{{read}}_fastq_screen.txt",
        png=f"{OUTDIR}/fastq_screen/fastq_screen_files/{{sample}}_{{lane}}_fq{{read}}_fastq_screen.png"
    threads: 
        get_resource("fastq_screen","threads")
    resources:
        mem_mb=get_resource("fastq_screen","mem_mb"),
        runtime=get_resource("fastq_screen","runtime")
    params:
        fastq_screen_config="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["parameters"]["fastq_screen_indexes"]["outdir"]),
        subset=100000,
        aligner='bowtie2'
    log:
        f"{LOGDIR}/fastq_screen_files/{{sample}}_{{lane}}_fq{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastq_screen_files/{{sample}}_{{lane}}_fq{{read}}.bmk"
    wrapper:
        "v1.23.4/bio/fastq_screen"


rule multiqc_files:
    input:
        multiqc_files_input(units)
    output:
        report(f"{OUTDIR}/multiqc/multiqc_files_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"]
    benchmark:
        f"{LOGDIR}/multiqc_files.bmk"
    log:
        f"{LOGDIR}/multiqc_files.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        runtime=get_resource("multiqc","runtime")
    conda:
        '../envs/multiqc.yaml'
    script:
        "../scripts/multiqc.py"


# QC for concat files

rule fastqc_concat:
    input:
        f"{OUTDIR}/trimmed/{{sample}}/{{sample}}_R{{read}}.fastq.gz"
    output:
        html=f"{OUTDIR}/qc/fastqc_concat/{{sample}}_R{{read}}_fastqc.html",
        zip_qc=f"{OUTDIR}/qc/fastqc_concat/{{sample}}_R{{read}}_fastqc.zip"
    threads: 
        get_resource("fastqc","threads")
    resources:
        mem_mb=get_resource("fastqc","mem_mb"),
        runtime=get_resource("fastqc","runtime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        f"{LOGDIR}/fastqc_concat/{{sample}}_R{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastqc_concat/{{sample}}_R{{read}}.bmk"
    wrapper:
        "0.74.0/bio/fastqc"


## Fastq_screen for concat files

rule fastq_screen_concat:
    input:
        f"{OUTDIR}/trimmed/{{sample}}/{{sample}}_R{{read}}.fastq.gz"
    output:
        txt=f"{OUTDIR}/fastq_screen/fastq_screen_concat/{{sample}}_R{{read}}_fastq_screen.txt",
        png=f"{OUTDIR}/fastq_screen/fastq_screen_concat/{{sample}}_R{{read}}_fastq_screen.png"
    threads: 
        get_resource("fastq_screen","threads")
    resources:
        mem_mb=get_resource("fastq_screen","mem_mb"),
        runtime=get_resource("fastq_screen","runtime")
    params:
        fastq_screen_config="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["parameters"]["fastq_screen_indexes"]["outdir"]),
        subset=100000,
        aligner='bowtie2'
    log:
        f"{LOGDIR}/fastq_screen_concat/{{sample}}_R{{read}}.log"
    benchmark:
        f"{LOGDIR}/fastq_screen_concat/{{sample}}_R{{read}}.bmk"
    wrapper:
        "v1.23.4/bio/fastq_screen"


rule multiqc_concat:
    input:
        multiqc_concat_input(units, final_step)
    output:
        report(f"{OUTDIR}/multiqc/multiqc_run_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"] + " --ignore fastqc_files"
    benchmark:
        f"{LOGDIR}/multiqc.bmk"
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        runtime=get_resource("multiqc","runtime")
    conda:
        '../envs/multiqc.yaml'
    script:
        "../scripts/multiqc.py"


# MultiQC for all executions

rule multiqc_join:
    input:
        OUTDIR
    output:
        report(f"{OUTDIR}/multiqc/multiqc_all_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"] + " --ignore fastqc_files"
    benchmark:
        f"{LOGDIR}/multiqc_all.bmk"
    log:
        f"{LOGDIR}/multiqc_all.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        runtime=get_resource("multiqc","runtime")
    conda:
        '../envs/multiqc.yaml'
    script:
        "../scripts/multiqc.py"
