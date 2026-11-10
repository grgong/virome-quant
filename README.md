# virome-quant

Counts the viruses of RNA-seq libraries by read mapping. The virus discovery (host removal,
assembly, taxonomy, vOTU clustering, CheckV) is done separately, from 2 M read pairs per library;
this repository makes that subsample and, with the resulting vOTU panel, counts the viruses in
all reads of every library. Needs SLURM.

## Setup

```bash
conda env create -p /scratch/$USER/envs/virome_quant -f environment.yml
cp config.sh my_config.sh        # set CONDA_ENV, CONDA_INIT, PARTITION, TMPDIR_JOBS
```

## Samplesheet

Tab-separated, one row per library, columns `sample`, `fastq_1`, `fastq_2`, `flowcell_lane`
(plus any columns of your own, copied to the results). Several lanes of one library:
comma-separated files in `fastq_1` / `fastq_2` and `FLOWCELL:LANE` items in `flowcell_lane`
(fields 3 and 4 of a read name: `zcat x_R1.fastq.gz | head -1 | cut -d: -f3,4`).
See `samples.example.tsv`.

## 1. Subsample (input of the discovery)

```bash
bash run.sh subsample my_config.sh samples.tsv /scratch/$USER/virome_2M
```

Draws 2,000,000 read pairs per library (seqtk, seed 11). Send the output directory
(`reads/`, `samples_2M.tsv`) to the discovery side.

## 2. Quantify

The discovery side returns a panel directory with `panel.fa`, `virus_units.tsv` and
`quality_summary.tsv`.

```bash
bash run.sh quantify my_config.sh samples.tsv <panel_dir> /scratch/$USER/virome_quant
```

All reads of every library -> fastp -> minimap2 to the panel -> CoverM. A vOTU counts in a
library when reads at >= 90 % identity cover >= 75 % of it; RPM = reads per million reads
passing fastp.

Results in `<outdir>/results/`: `viral_rpm.tsv` (libraries x high-quality vOTUs, RPM),
`votus.tsv` (every vOTU), `viral_long.tsv` (library x vOTU: reads, covered fraction, RPM, carry-over
columns).

A failed job: read its log in `<outdir>/logs/`, fix the cause, then `bash run.sh resume <outdir>`.
