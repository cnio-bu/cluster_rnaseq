import pandas as pd
import numpy as np

samples = pd.read_table(snakemake.params.samples).set_index("sample", drop=False)

counts = [pd.read_table(f, index_col=0, usecols=[0, 6], header=None, skiprows=2)
          for f in snakemake.input]

for t, sample in zip(counts, samples.index):
    t.columns = [sample]

count_matrix = pd.concat(counts, axis=1)
count_matrix.index.name = "geneID"

# collapse technical replicates
count_matrix = count_matrix.groupby(count_matrix.columns, axis=1).sum()
count_matrix = count_matrix.apply(np.floor)
count_matrix.to_csv(snakemake.output[0], sep="\t")