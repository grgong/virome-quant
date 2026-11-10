#!/usr/bin/env bash
#SBATCH --job-name=vq_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# Map all reads of each library to the vOTU panel, one array task per library:
#   sbatch --array=1-N map.sbatch samples.tsv panel_dir outdir
set -euo pipefail
SAMPLES=$1; PANEL=$2/panel.fa; OUT=$3

IFS=$'\t' read -r sample fq1 fq2 < <(awk -F'\t' -v OFS='\t' -v n="$SLURM_ARRAY_TASK_ID" '
    NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    NR == n + 1 {print $c["sample"], $c["fastq_1"], $c["fastq_2"]}' "$SAMPLES")

mkdir -p "$OUT/map" "$OUT/fastp"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# fastp needs real files, so join the lanes first
if [[ $fq1 == *,* ]]; then
    cat ${fq1//,/ } > "$tmp/R1.fq.gz"; cat ${fq2//,/ } > "$tmp/R2.fq.gz"
    fq1=$tmp/R1.fq.gz; fq2=$tmp/R2.fq.gz
fi

fastp -i "$fq1" -I "$fq2" -o "$tmp/T1.fq.gz" -O "$tmp/T2.fq.gz" -z 1 \
    --detect_adapter_for_pe --trim_poly_g --low_complexity_filter --length_required 50 \
    -w "$SLURM_CPUS_PER_TASK" -j "$OUT/fastp/$sample.json" -h /dev/null 2> /dev/null

# keep primary alignments that cover at least 80 % of the read
minimap2 -ax sr -t "$SLURM_CPUS_PER_TASK" "$PANEL" "$tmp/T1.fq.gz" "$tmp/T2.fq.gz" 2> /dev/null |
    samtools view -u -F 0x904 -e '(qlen-sclen) >= 0.8*qlen' - |
    samtools sort -@ 2 -o "$tmp/out.bam" -
mv "$tmp/out.bam" "$OUT/map/$sample.bam"
