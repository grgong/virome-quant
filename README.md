# virome-quant

Count viruses in RNA-seq libraries by mapping the reads to a set of viral genomes (vOTUs) that
were assembled separately. Needs SLURM.

The virus discovery (host removal, assembly, taxonomy, vOTU clustering, CheckV) runs on 2 M read
pairs per library and is done elsewhere. This repository does the two steps around it: it makes
the 2 M subsample, and once you have the vOTU panel it counts the vOTUs in all reads.

## Setup

```bash
git clone https://github.com/grgong/virome-quant.git
cd virome-quant
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
N=$(( $(wc -l < samples.tsv) - 1 ))      # number of samples
sbatch -p <partition> --array=1-$N%20 subsample.sbatch samples.tsv virome_2M
```

Send `virome_2M/reads/` (about 280 MB per library) to Gaorui, who runs the discovery.

## Count

Gaorui sends back a directory with `panel.fa` (the vOTU sequences), `virus_units.tsv` and
`quality_summary.tsv`.

Counting uses all reads of each library, not the 2 M subsample. `map.sbatch` trims the reads with
fastp and maps them to the vOTU sequences with minimap2, keeping the reads aligned over at least
80 % of their length. `count.py` then runs CoverM on the BAM files. A read counts for a vOTU when
it is at least 90 % identical to it. A vOTU is detected in a library when such reads cover at
least 75 % of its length; below that its count is set to 0. RPM is the count per million reads
that pass fastp.

```bash
conda activate virome_quant
export TMPDIR=/scratch/$USER/tmp
N=$(( $(wc -l < samples.tsv) - 1 ))      # number of samples, as for the subsample
sbatch -p <partition> --array=1-$N%20 map.sbatch samples.tsv <panel_dir> virome_quant
python3 count.py samples.tsv <panel_dir> virome_quant      # after every map task has finished
```

If a map task fails, rerun it alone with `--array=<its number>`.

Results are in `virome_quant/results/`:

- `viral_long.tsv`: one row per library and vOTU, with reads, fraction covered, RPM and CheckV quality
- `viral_rpm.tsv`: libraries by high-quality vOTUs, RPM

## Plot

```bash
Rscript plot.R virome_quant/results/viral_rpm.tsv -s samples_to_plot.txt -g condition -o figure
```

Draws a heatmap with the vOTUs in rows, the samples in columns and the RPM in five bins, to
`figure.pdf` and `figure.png`. It needs R and ComplexHeatmap, both in the environment.

- `-s`: a text file with one sample name per line. Only these samples are drawn, in this order,
  and vOTUs without reads in them are left out. Without `-s` every sample is drawn.
- `-g`: a column of the samples table, such as `condition`. The samples are coloured by it, and
  neighbouring samples with the same value form a group.
- `-o`: output prefix.

Cells are 3.2 mm squares. The script says so when the figure is wider than 183 mm.
