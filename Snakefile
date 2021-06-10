import glob
import os
import pandas as pd
import numpy as np
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


def get_aligner(chosen_aligner:int) -> str:

	available_aligners = {0:'star', 1:'salmon', 2:'hisat2'}
	try:
		return available_aligners[chosen_aligner]
	except KeyError:
		print(f'Invalid aligner choice: {chosen_aligner}')
		sys.exit(1)
	

#### LOAD SAMPLES TABLES ###

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")


#### Load rules ####
include: 'rules/common.smk'
include: 'rules/qc.smk'
include: 'rules/preprocess.smk'
include: 'rules/index.smk'
include: 'rules/align.smk'
include: 'rules/quantification.smk'

rule all:
	input:
		f"{OUTDIR}/qc/multiqc_report.html",
		f"{OUTDIR}/qc_concat/multiqc_report.html",
		expand(f'{OUTDIR}/quant/{chosen_aligner}/{{sample}}.tab', sample=samples['sample'])
