# virome-quant

Count viruses in RNA-seq libraries by mapping the reads to a set of viral genomes (vOTUs) that
were assembled separately. Needs SLURM.

The virus discovery (host removal, assembly, taxonomy, vOTU clustering, CheckV) runs on 2 M read
pairs per library and is done elsewhere. This repository does the two steps around it: it makes
the 2 M subsample, and once you have the vOTU panel it counts the vOTUs in all reads.

## Setup

```bash
conda env create -p ~/envs/virome_quant -f environment.yml
conda activate virome_quant
export TMPDIR=/scratch/$USER/tmp; mkdir -p $TMPDIR    # a mapping job uses up to 20 GB here
```

Activate the environment in the shell you submit from; the jobs inherit it.

## Samplesheet

One row per sample, tab-separated with Unix line endings, columns `sample`, `fastq_1` (R1 file) and
`fastq_2` (R2 file), with absolute paths. Extra columns, such as condition, are copied to the
results. See `samples.example.tsv`.

## Subsample

```bash
N=$(( $(wc -l < samples.tsv) - 1 ))
sbatch -p <partition> --array=1-$N%20 subsample.sbatch samples.tsv virome_2M
```

Send `virome_2M/reads/` (about 280 MB per library) to whoever runs the discovery.

## Count

The discovery gives you a directory with `panel.fa`, `virus_units.tsv` and `quality_summary.tsv`.

```bash
sbatch -p <partition> --array=1-$N%20 map.sbatch samples.tsv <panel_dir> virome_quant
python3 count.py samples.tsv <panel_dir> virome_quant      # after every map task has finished
```

A vOTU counts in a library when reads at least 90 % identical to it cover at least 75 % of its
length. RPM is per million reads that pass fastp. If a map task fails, rerun it alone with
`--array=<its number>`.

Results are in `virome_quant/results/`:

- `viral_long.tsv`: one row per library and vOTU, with reads, fraction covered, RPM and CheckV quality
- `viral_rpm.tsv`: libraries by high-quality vOTUs, RPM

## Plot

```bash
python3 plot.py virome_quant/results/viral_rpm.tsv -s samples_to_plot.txt -g condition -o figure
```

Draws a heatmap with the vOTUs in rows, the samples in columns and the RPM in five bins, to
`figure.pdf` and `figure.png`.

- `-s`: a text file with one sample name per line. Only these samples are drawn, in this order,
  and vOTUs without reads in them are left out. Without `-s` every sample is drawn.
- `-g`: a column of the samples table, such as `condition`. The samples are coloured by it, and
  neighbouring samples with the same value form a group.
- `-o`: output prefix.

Cells are 3.2 mm squares. The script says so when the figure is wider than 183 mm.
