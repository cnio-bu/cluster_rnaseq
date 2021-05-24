import hashlib
import os

import pandas as pd
import numpy  as np

def get_md5(read):

    cur_dir   = os.getcwd()
    read_path = os.path.join(cur_dir, read)

    if os.path.isfile(read_path):
        with open(read_path, "rb") as file_name:
                content = file_name.read()
                aux_var = hashlib.md5()
                aux_var.update(content)
        return aux_var.hexdigest()
    else:
        return 'missing_file'


## SNAKEMAKE I/O ##
units = snakemake.input['units']

## SNAKEMAKE PARAMS ##
logdir = snakemake.params['logdir']


units  = pd.read_csv(units, sep='\t')

md5s    = units.set_index('fq1')['md5_fq1'].to_dict()
md5s.update(units.set_index('fq2')['md5_fq2'].to_dict())

md5_check = pd.DataFrame.from_dict(md5s, orient='index', columns=['md5'])
md5_check['fastq'] = md5_check.index
md5_check = md5_check.reindex(columns=['fastq', 'md5'])

md5_check = md5_check.dropna(subset=["fastq", "md5"])

md5_check['md5_fq_check'] = md5_check[md5_check['fastq'].notnull()]['fastq'].apply(get_md5)

if not md5_check.empty:
    md5_check['md5_fq_check'] = np.where(md5_check['md5_fq_check'].str.lower() == md5_check['md5'].str.lower(), 'Passed', 'Failed')
    md5_check[md5_check['md5'].notnull()].to_csv(f"{logdir}/md5.done.tsv", sep = '\t', index=False)

if (md5_check[md5_check['md5'].notnull()]['md5_fq_check'] == 'Failed').any():
    sys.exit(f"\nSome files did not pass the checksum veryfication.\nCheck the report at {logdir}/md5.done.tsv\nExiting...\n")
else:
    sys.exit(0)
