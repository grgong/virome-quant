# virome-quant

Scripts to measure how much of each virus there is in RNA-seq libraries. The reads are mapped to
a set of viral genomes (vOTUs) and counted. You need a cluster with SLURM.

The viral genomes are not built here. I assemble them from 2 M read pairs of each library
(host removal, assembly, taxonomy, vOTU clustering, CheckV). The scripts here do what comes
before and after: first they draw the 2 M read pairs, and when the viral genomes are ready they
count them in all the reads.

## Setup

```bash
git clone https://github.com/grgong/virome-quant.git
cd virome-quant
conda env create -p ~/envs/virome_quant -f environment.yml
conda activate virome_quant
export TMPDIR=/scratch/$USER/tmp; mkdir -p $TMPDIR    # the mapping jobs sort their BAM files here
```

Submit the jobs from a shell where the environment is active. They use the tools of that shell.

## Samplesheet

A tab-separated file with one line per sample and three columns: `sample`, `fastq_1` (the R1
file) and `fastq_2` (the R2 file). Give the full paths of the files. You can add other columns,
for example `condition`; they are copied to the result tables. `samples.example.tsv` is an
example.

If you make the file on Windows, run `dos2unix samples.tsv` on the cluster before using it.

## Subsample

```bash
N=$(( $(wc -l < samples.tsv) - 1 ))      # number of samples
sbatch -p <partition> --array=1-$N%20 subsample.sbatch samples.tsv virome_2M
```

This writes 2 M read pairs of each library to `virome_2M/reads/`, about 280 MB per library. Send
that directory to me.

## Count

I send back a directory with three files: `panel.fa` (the viral genomes),
`virus_units.tsv` and `quality_summary.tsv`.

The counting uses all the reads of each library, not the 2 M subsample. `map.sbatch` trims the
reads with fastp and maps them to the viral genomes with minimap2. A read is kept when at least
80 % of its length is aligned. `count.py` then runs CoverM on the BAM files and counts, for each
vOTU, the reads that are at least 90 % identical to it. A vOTU is called detected in a library
when these reads cover at least 75 % of its length. If they cover less, its count is set to 0.
RPM is the count per million reads that pass fastp.

```bash
conda activate virome_quant
export TMPDIR=/scratch/$USER/tmp
N=$(( $(wc -l < samples.tsv) - 1 ))      # number of samples, the same as for the subsample
sbatch -p <partition> --array=1-$N%20 map.sbatch samples.tsv <panel_dir> virome_quant
python3 count.py samples.tsv <panel_dir> virome_quant      # when all the map tasks have finished
```

If one map task fails, submit it again alone with `--array=<its number>`.

The results are in `virome_quant/results/`:

- `viral_long.tsv` has one line per library and vOTU, with the number of reads, the fraction of
  the genome covered, the RPM and the CheckV quality.
- `viral_rpm.tsv` has the RPM with one line per library and one column per high-quality vOTU.

## Plot

```bash
Rscript plot.R virome_quant/results/viral_rpm.tsv -s samples_to_plot.txt -g condition -o figure
```

This draws a heatmap of `viral_rpm.tsv` and saves it as `figure.pdf` and `figure.png`. The vOTUs
are the rows, the samples are the columns, and the colour gives the RPM in five classes. R and
ComplexHeatmap are in the environment.

- `-s` is a text file with one sample name per line. Only these samples are drawn, in the order
  of the file. A vOTU with no reads in these samples is not drawn. Without `-s`, all the samples
  are drawn.
- `-g` is a column of the samplesheet, for example `condition`. The samples get the colour of
  their value in this column, and samples next to each other with the same value are drawn as
  one group.
- `-o` is the name of the output files, without the extension.
