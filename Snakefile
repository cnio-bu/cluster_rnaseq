import glob
import os
import re
import pandas as pd
import numpy as np
import sys
from snakemake.utils import validate, min_version

#### GLOBAL PARAMETERS ####

min_version('6.2.1')

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
IS_GENCODE = ('--gencode' if config['parameters']['salmon_index']['gencode'] else '')

#### GLOBAL scope functions ####
def get_resource(rule,resource) -> int:
	'''
	Attempt to parse config.yaml to retrieve resources available for a given
	rule. It will revert to default if a key error is found. Returns an int.
	with the allocated resources available for said rule. Ex: "threads": 1
	'''

	try:
		return config['resources'][rule][resource]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]

#TODO: add defaults
def get_params(rule,param) -> int:
	'''
	Attempt to parse config.yaml to retrieve parameters available for a given
	rule. It will crash otherwise.
	''' 
	try:
		return config['parameters'][rule][param]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve parameter for {rule}/{param}: Exiting...')
		sys.exit(1)


def get_aligner(chosen_aligner:int, UMIs) -> str:

	if 	UMIs: 
		available_aligners = {0:'star', 2:'hisat2'}
	else:
		available_aligners = {0:'star', 1:'salmon', 2:'hisat2'}
	try:
		return available_aligners[chosen_aligner]
	except KeyError:
		print(f'Invalid aligner choice: {chosen_aligner}. Remember if ' +
		'your RNAseq contains UMIs, only star(0) or hisat2(2) could be chosen)')
		sys.exit(1)


def get_quantifier(chosen_quantifier:int) -> str:

	available_quantifiers = {0:'htseq', 1:'featureCounts'}
	try:
		return available_quantifiers[chosen_quantifier]
	except KeyError:
		print(f'Invalid quantifier choice: {chosen_quantifier}')
		sys.exit(1)

	
def get_reference_level(x):
    '''
	Gets the reference level for a column.
	'''
    levels = x.unique()
    is_reference = [lvl.startswith("*") for lvl in levels]
    if any(is_reference):    				
    	if sum(is_reference) > 1:
    		raise ValueError(f'More than two reference levels were specified ' +
							 f'for the same column: {levels[is_reference]}')
    	reference = levels[is_reference][0].strip("*")
    else:
    	reference = sorted(levels)[0]
    return reference


def get_covariates(design):
	'''
	Get the covariates and the variable of interest (the last one) for differential expression.
	'''
	covariates = list(filter(None, re.split("[ \+\*:~]", design)))
	variable_interest = covariates[-1]
	covariates = list(set(covariates))
	return {"covariates": covariates, "variable_interest": variable_interest}


def get_final_step():
	arg_list = sys.argv
	possible_steps = ["all", "index", "files_qc", "trimming", "alignment", "quantification", \
					"diffexp", "plots", "multiqc_all"]

	step = [i for i in arg_list if i in possible_steps]

	if len(step) > 1:
		print(f'Failed to resolve target rule for {step}. Only one can be specified. Exiting...')
		sys.exit(1)

	elif len(step) == 0:
		step = ["plots"]
		
	elif step[0] == "all":
		step[0] = "plots"

	return step[0]



#### LOAD SAMPLES TABLES ###

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

designmatrix = pd.read_table(config["parameters"]["deseq2"]["designmatrix"], \
							 dtype=str).set_index("sample", drop=False)
validate(designmatrix, schema="schemas/designmatrix.schema.yaml")

#### Get information about UMI analysis.
UMIs = get_params("umi_processing", "enabled")

#### Get aligner ####
chosen_aligner = get_aligner(int(config['aligner']), UMIs)

#### Get quantifier ####
chosen_quantifier = get_quantifier(int(config['quantifier']))

#### Auxiliar variable to generate the paths after quantification, due to salmon does not need to 
#### quantify after align while STAR and hisat2 do, so they need an extra folder in the results.
deseq_path = f"{chosen_aligner}" if chosen_aligner == "salmon" else f"{chosen_aligner}/{chosen_quantifier}"

#### Subset the design matrix keeping only samples for DEA ####
DEAsamples = samples[samples["diffexp"]].index
designmatrix = designmatrix.loc[DEAsamples,]

#### Get reference level for each column in design matrix ####
designmatrix.drop(columns="sample", inplace=True)
ref_levels = designmatrix.apply(get_reference_level, axis="index")

#### Remove '*' prefix from design matrix cells ####
designmatrix = designmatrix.apply(lambda row: [str(x).removeprefix("*") \
                                              for x in row])

#### Get column of interest for differential expression
var_info = get_covariates(config["parameters"]["deseq2"]["design"])
var_interest = var_info["variable_interest"]
covariates = var_info["covariates"]

#### Contrasts ####
ref_interest = ref_levels[var_interest]
rest_levels = [y for y in designmatrix[var_interest].unique() \
              if y != ref_interest]
contrasts = [(z + "_vs_" + ref_interest, [z, ref_interest]) \
            for z in rest_levels]
contrasts = {key: value for (key, value) in contrasts}

if len(rest_levels) > 1:
	allSamples = {"allSamples": list(set([ref_interest] + rest_levels))}
	allSamples.update(contrasts)
else:
	allSamples = contrasts

### Batch correction (for plotting PCAs and correlations)
filesuffix = [""]
if len(covariates) > 1:
	filesuffix += ["_batchCorrected"]
	batch = [x for x in covariates if x != var_interest]
else:
	batch = None


#### Final Step #### Global variable for the last step in the run
final_step = get_final_step()


#### Load rules ####
include: 'rules/common.smk'
include: 'rules/qc.smk'
include: 'rules/preprocess.smk'
include: 'rules/index.smk'
include: 'rules/align.smk'
include: 'rules/quantification.smk'
include: 'rules/deseq2.smk'
include: 'rules/plots.smk'


def get_index_input():

	index_input = config["ref"][chosen_aligner][f"{chosen_aligner}_index"]
	return index_input


def get_trimming_input():

	trimming_input = [f"{OUTDIR}/multiqc/multiqc_files_report.html"]
	
	for sample in samples['sample']:
		if single_end:
			trimming_input += [f"{OUTDIR}/qc/fastqc_concat/{sample}_R1_fastqc.html"]
		else:
			trimming_input += expand(f"{OUTDIR}/qc/fastqc_concat/{sample}_R{{strand}}_fastqc.html", strand=[1,2])

	## Deal with optional rules (fastq_screen). If fastq_screen is enabled, this tool will also be performed over 
	# concatenated files.
	if config["parameters"]["fastq_screen"]["enabled"]:
		for sample in samples['sample']:
			if single_end:
				trimming_input += [f"{OUTDIR}/fastq_screen/fastq_screen_concat/{sample}_R1_fastq_screen.txt"]
			else:
				trimming_input += expand(f"{OUTDIR}/fastq_screen/fastq_screen_concat/{sample}_R{{strand}}_fastq_screen.txt", strand=[1,2])

	for sample in samples['sample']:
		if single_end:
			trimming_input += [f"{OUTDIR}/trimmed/{sample}/{sample}_R1.fastq.gz"]
		else:
			trimming_input += expand(f"{OUTDIR}/trimmed/{sample}/{sample}_R{{strand}}.fastq.gz", strand=[1,2])
	
	#trimming_input += [f"{OUTDIR}/multiqc/multiqc_run_report.html"]
	return trimming_input


def get_alignment_input():

	alignment_input = [f"{OUTDIR}/multiqc/multiqc_files_report.html"]
	
	for sample in samples['sample']:
		if single_end:
			alignment_input += [f"{OUTDIR}/qc/fastqc_concat/{sample}_R1_fastqc.html"]
		else:
			alignment_input += expand(f"{OUTDIR}/qc/fastqc_concat/{sample}_R{{strand}}_fastqc.html", strand=[1,2])
	
	if config["parameters"]["fastq_screen"]["enabled"]:
		for sample in samples['sample']:
			if single_end:
				alignment_input += [f"{OUTDIR}/fastq_screen/fastq_screen_concat/{sample}_R1_fastq_screen.txt"]
			else:
				alignment_input += expand(f"{OUTDIR}/fastq_screen/fastq_screen_concat/{sample}_R{{strand}}_fastq_screen.txt", strand=[1,2])
	
	if chosen_aligner == "salmon":
		alignment_input += expand(f"{OUTDIR}/quant/salmon/{{sample}}/quant.sf", sample=samples['sample'])
	else:
		#if UMIs:
			#alignment_input += expand(f"{OUTDIR}/dedup/alignments/{{sample}}/deduplicated.bam", sample=samples['sample'])
		
		alignment_input += expand(f"{OUTDIR}/mapped/{chosen_aligner}/{{sample}}/Aligned.sortedByCoord.out.bam", sample=samples['sample'])
	#alignment_input += [f"{OUTDIR}/multiqc/multiqc_run_report.html"]
	return alignment_input


def get_quantification_input():

	quantification_input= [f"{OUTDIR}/multiqc/multiqc_files_report.html"]
	
	quantification_input += [f"{OUTDIR}/deseq2/{deseq_path}/counts.tsv"]
	
	quantification_input += [f"{OUTDIR}/multiqc/multiqc_run_report.html"]
	return quantification_input


def get_diffexp_input():

	diffexp_input = [f"{OUTDIR}/multiqc/multiqc_files_report.html"]

	diffexp_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.xlsx", contrast=contrasts.keys())
	diffexp_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.tsv", contrast=contrasts.keys())
	
	diffexp_input += [f"{OUTDIR}/multiqc/multiqc_run_report.html"]
	return diffexp_input


def get_plots_input():

	plots_input= [f"{OUTDIR}/multiqc/multiqc_files_report.html"]
	
	plots_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.xlsx", contrast=contrasts.keys())
	plots_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/{{contrast}}_diffexp.tsv", contrast=contrasts.keys())
	plots_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_topbottomDEgenes.{{pext}}", \
						contrast=contrasts.keys(), pext = ["pdf", "png"])
	plots_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{ALLcontrast}}/plots/{{ALLcontrast}}_{{plot}}{{fsuffix}}.{{pext}}", \
						ALLcontrast=allSamples.keys(), fsuffix=filesuffix, plot = ["pca", "dist"], pext = ["pdf", "png"])
	plots_input += expand(f"{OUTDIR}/deseq2/{deseq_path}/{{contrast}}/plots/{{contrast}}_MAplot.{{pext}}", \
						contrast=contrasts.keys(), pext = ["pdf", "png"])

	plots_input += [f"{OUTDIR}/multiqc/multiqc_run_report.html"]
	return plots_input



# TARGET RULES
rule all:
	input:
		get_plots_input()


rule index:
	input:
		get_index_input()


# Check files' MD5 and creates a MultiQC report using the fastqc's reports for original files
rule files_qc:
	input:
		f"{OUTDIR}/multiqc/multiqc_files_report.html"


rule trimming:
	input:
		get_trimming_input()


rule alignment:
	input:
		get_alignment_input()


rule quantification:
	input:
		get_quantification_input()


rule diffexp:
	input:
		get_diffexp_input()


rule plots:
	input:
		get_plots_input()


# Creates a MultiQC report for all files that it founds, mixing all aligners or quantifiers that has been ran
rule multiqc_all:
	input:
		f"{OUTDIR}/multiqc/multiqc_all_report.html"
