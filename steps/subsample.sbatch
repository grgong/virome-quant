#!/usr/bin/env bash
# Subsample, one array task per sample: 2,000,000 random read pairs (seqtk, seed 11). These are
# the reads the virus discovery (host removal, assembly, classification) uses; it is run on the
# discovery side with the same draw, so it gets the same reads as from the full files. Lanes are
# read in samplesheet order; R1 and R2 are drawn with the same seed, which keeps the pairs in sync
# (seqtk reservoir sampling depends only on the seed and the number of records).
# Outputs: reads/<sample>_R{1,2}.fq.gz, reads/<sample>.md5
source "${VP}/steps/common.sh"
N_PAIRS=2000000
SEED=11

IFS=$'\t' read -r sample fq1 fq2 <<< "$(sample_row "${SLURM_ARRAY_TASK_ID}")"
[ -e "${OUT}/done/subsample.${sample}" ] && { log "${sample} already done"; exit 0; }
mkdir -p "${OUT}/reads"
log "${sample}"

cat ${fq1//,/ } | seqtk sample -s"${SEED}" - "${N_PAIRS}" > "${SCR}/R1.fq"
cat ${fq2//,/ } | seqtk sample -s"${SEED}" - "${N_PAIRS}" > "${SCR}/R2.fq"
name() { awk -v line="$2" 'NR==line {sub(/^@/, ""); sub(/[ \/].*/, ""); print; exit}' "$1"; }
n1=$(( $(wc -l < "${SCR}/R1.fq") / 4 )); n2=$(( $(wc -l < "${SCR}/R2.fq") / 4 ))
[ "${n1}" -gt 0 ] && [ "${n1}" -eq "${n2}" ] || { echo "ERROR: R1/R2 record counts ${n1} / ${n2}" >&2; exit 1; }
last=$(( (n1 - 1) * 4 + 1 ))
[ "$(name "${SCR}/R1.fq" 1)" = "$(name "${SCR}/R2.fq" 1)" ] &&
[ "$(name "${SCR}/R1.fq" "${last}")" = "$(name "${SCR}/R2.fq" "${last}")" ] ||
    { echo "ERROR: R1/R2 out of sync after subsampling" >&2; exit 1; }

for r in 1 2; do
    pigz -p "${THREADS}" -c "${SCR}/R${r}.fq" > "${SCR}/R${r}.fq.gz"
    mv "${SCR}/R${r}.fq.gz" "${OUT}/reads/${sample}_R${r}.fq.gz"
done
( cd "${OUT}/reads" && md5sum "${sample}_R1.fq.gz" "${sample}_R2.fq.gz" > "${sample}.md5" )
touch "${OUT}/done/subsample.${sample}"
log "${sample}: ${n1} pairs"
