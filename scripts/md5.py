import hashlib
import pandas as pd
import numpy as np
import os

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


# Main script
units = snakemake.params.units_df
LOGDIR = snakemake.params.log_path

units_2 = units.copy()

md5s    = units_2.set_index('fq1')['md5_fq1'].to_dict()
md5s.update(units_2.set_index('fq2')['md5_fq2'].to_dict())

md5_check = pd.DataFrame.from_dict(md5s, orient='index', columns=['md5'])
md5_check['fastq'] = md5_check.index
md5_check = md5_check.reindex(columns=['fastq', 'md5'])

md5_check['md5_fq_check'] = md5_check[md5_check['fastq'].notnull()]['fastq'].apply(get_md5)

md5_check['md5_fq_check'] = np.where(md5_check['md5_fq_check'] == md5_check['md5'], 'Passed', 'Failed')

#print(md5_check)
md5_check[md5_check['md5'].notnull()].to_csv(f"{LOGDIR}/md5.done.tsv", sep = '\t', index=False)

if (md5_check[md5_check['md5'].notnull()]['md5_fq_check'] == 'Failed').any():
    sys.exit(f"\nSome files did not pass the checksum veryfication.\nCheck the report at {LOGDIR}/md5.done.tsv\nExiting...\n")
else:
    sys.exit(0)
