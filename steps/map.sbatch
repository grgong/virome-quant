#!/usr/bin/env bash
# Mapping, one array task per sample: all reads -> fastp (same options as on the discovery side)
# -> minimap2 -x sr to the vOTU panel; primary alignments with >= 80 % of the trimmed read
# aligned are kept. No host removal here. Reads total = reads passing fastp (both mates), from
# its JSON. fastp needs real files (from a pipe it reads 0 reads and exits 0), so the lanes of a
# multi-lane sample are concatenated to scratch first.
# Outputs: map/<sample>.bam(.bai), map/<sample>.nreads, fastp/<sample>.json
source "${VP}/steps/common.sh"
PANEL="${OUT}/panel/panel.fa"
FASTP_ARGS=(--detect_adapter_for_pe --trim_poly_g --low_complexity_filter --length_required 50)

IFS=$'\t' read -r sample fq1 fq2 <<< "$(sample_row "${SLURM_ARRAY_TASK_ID}")"
[ -e "${OUT}/done/map.${sample}" ] && { log "${sample} already done"; exit 0; }
mkdir -p "${OUT}/map" "${OUT}/fastp"
log "${sample}"

if [[ "${fq1}" == *,* ]]; then   # lanes in samplesheet order, same for R1 and R2 (concatenated gzip is valid gzip)
    cat ${fq1//,/ } > "${SCR}/R1.fq.gz"; cat ${fq2//,/ } > "${SCR}/R2.fq.gz"
    in1="${SCR}/R1.fq.gz"; in2="${SCR}/R2.fq.gz"
else
    in1="${fq1}"; in2="${fq2}"
fi
# --trim_poly_g explicit: fastp 1.0.1 detects NovaSeq X reads as two-colour only from 1.1.0 on
fastp -i "${in1}" -I "${in2}" -o "${SCR}/T1.fq.gz" -O "${SCR}/T2.fq.gz" -z 1 "${FASTP_ARGS[@]}" \
    -w "${THREADS}" -j "${OUT}/fastp/${sample}.json" -h "${SCR}/fastp.html" 2> "${SCR}/fastp.log"
rm -f "${SCR}/R1.fq.gz" "${SCR}/R2.fq.gz"
read -r raw passed < <(python3 -I -c 'import json, sys; s = json.load(open(sys.argv[1]))["summary"]; print(s["before_filtering"]["total_reads"], s["after_filtering"]["total_reads"])' \
    "${OUT}/fastp/${sample}.json")
[ "${raw}" -gt 0 ] && [ "${passed}" -gt 0 ] || { echo "ERROR: fastp read ${raw} reads, passed ${passed}" >&2; exit 1; }

minimap2 -ax sr -t "${THREADS}" "${PANEL}" "${SCR}/T1.fq.gz" "${SCR}/T2.fq.gz" 2> "${SCR}/mm2.log" |
    samtools view -u -F 0x904 -e '(qlen-sclen) >= 0.8*qlen' - |
    samtools sort -@ 2 -o "${SCR}/out.bam" -
mapped=$(grep -oE 'mapped [0-9]+ sequences' "${SCR}/mm2.log" | awk '{s += $2} END {print s}')
[ "${mapped}" -eq "${passed}" ] || { echo "ERROR: minimap2 read ${mapped} reads, fastp passed ${passed}" >&2; exit 1; }
mv "${SCR}/out.bam" "${OUT}/map/${sample}.bam"
samtools index "${OUT}/map/${sample}.bam"
printf 'sample\treads_total\treads_raw\n%s\t%s\t%s\n' "${sample}" "${passed}" "${raw}" > "${OUT}/map/${sample}.nreads"
touch "${OUT}/done/map.${sample}"
log "${sample}: ${raw} raw reads, ${passed} after fastp, $(samtools view -c "${OUT}/map/${sample}.bam") on the panel"
