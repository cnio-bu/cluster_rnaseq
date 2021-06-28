# This is the same code that STAR uses in its wrapper, but we change zcat for gunzip due to some problems
# running the pipeline in mac's zsh terminal

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])

if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand gunzip -c"
else:
    readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"

shell(
    "STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.input.index} "
    "--readFilesIn {input_str} "
    "{readcmd} "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "{log}"
)
