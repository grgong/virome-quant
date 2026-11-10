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

Tab-separated, Unix line endings, absolute paths, columns `sample`, `fastq_1`, `fastq_2`. Other
columns (condition, stage, ...) are copied to the results. A library on several lanes has the
files separated by commas in `fastq_1` and `fastq_2`, in the same order. See `samples.example.tsv`.

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
