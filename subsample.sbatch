#!/usr/bin/env bash
#SBATCH --job-name=vq_subsample
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# Take 2 M read pairs from each library, one array task per library:
#   sbatch --array=1-N subsample.sbatch samples.tsv outdir
set -euo pipefail
SAMPLES=$1; OUT=$2

IFS=$'\t' read -r sample fq1 fq2 < <(awk -F'\t' -v OFS='\t' -v n="$SLURM_ARRAY_TASK_ID" '
    NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    NR == n + 1 {print $c["sample"], $c["fastq_1"], $c["fastq_2"]}' "$SAMPLES")

mkdir -p "$OUT/reads"
# same seed for R1 and R2, so the pairs stay in sync
seqtk sample -s11 "$fq1" 2000000 | gzip > "$OUT/reads/${sample}_R1.fq.gz"
seqtk sample -s11 "$fq2" 2000000 | gzip > "$OUT/reads/${sample}_R2.fq.gz"
cd "$OUT/reads" && md5sum "${sample}"_R?.fq.gz > "$sample.md5"
