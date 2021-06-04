 ## SALMON RULES ##

rule generate_decoy_sequences:
    input:
        genome=config['ref']['salmon']['genome_assembly']
    output:
        decoys="/".join(config["ref"]["salmon"]["salmon_index"].split('/')[:-1])+"/decoys.txt"
    resources:
        walltime=1
    shell: 
        """
        grep "^>" <(gunzip -c {input.genome}) | cut -d " "  -f 1 > {output.decoys}
        sed -i -e 's/>//g' {output.decoys}
        """


rule build_gentrome:
    input:
        genome=config['ref']['salmon']['genome_assembly'],
        transcriptome=config['ref']['salmon']['transcriptome']
    output:
        gentrome="/".join(config["ref"]["salmon"]["salmon_index"].split('/')[:-1])+"/gentrome.fa.gz"
    resources:
        walltime=1
    shell:
        'cat {input.transcriptome} {input.genome} > {output.gentrome}'


rule salmon_index:
    input:
        gentrome=rules.build_gentrome.output.gentrome,
        decoys=rules.generate_decoy_sequences.output.decoys
    output:
        directory(config["ref"]["salmon"]["salmon_index"])
    threads:
        get_resource('salmon_index', 'threads')
    resources:
        mem=get_resource('salmon_index', 'mem'),
        walltime=get_resource('salmon_index', 'walltime')
    params:
        gencode = IS_GENCODE # Dirty, but necessary. See Snakefile.
    conda:
        '../envs/aligners.yaml'
    shell:
        'salmon index -t {input.gentrome} -d {input.decoys} -p {threads} -i {output} {params.gencode}'


## STAR RULES
rule star_index:
    input:
        fasta = config["ref"]["star"]["fasta"] if config["ref"]["star"]["fasta"] else "-",
        gtf   = config["ref"]["star"]["annotation"] if config["ref"]["star"]["annotation"] else "-"
    output:
        directory(config["ref"]["star"]["star_index"])
    threads: 
        get_resource("star_index", "threads")
    resources:
        mem=get_resource("star_index", "mem"),
        walltime=get_resource("star_index", "walltime")
    params:
        extra = ""
    log:
        f"{LOGDIR}/star_index/index.log"
    wrapper:
        "0.74.0/bio/star/index"


## HISAT-2 RULES
rule hisat2_index:
    input:
        fasta = config["ref"]["hisat2"]["fasta"] if config["ref"]["hisat2"]["fasta"] else "-"
    output:
        directory(config["ref"]["hisat2"]["hisat2_index"])
    threads:
        get_resource("hisat2_index", "threads")
    resources:
        mem=get_resource("hisat2_index", "mem"),
        walltime=get_resource("hisat2_index", "walltime")
    params:
        prefix=config["ref"]["hisat2"]["hisat2_index"] + "/hisat2_index"
    log:
        f"{LOGDIR}/hisat2_index/index.log"
    wrapper:
        "0.74.0/bio/hisat2/index"
