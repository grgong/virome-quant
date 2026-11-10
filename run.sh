#!/usr/bin/env bash
# Virome quantification: checks the inputs, then submits the steps as SLURM jobs. Run on the login node.
#
#   bash run.sh subsample <config.sh> <samples.tsv> <outdir>
#        2 M read pairs per sample, the input of the virus discovery (sent to the discovery side)
#   bash run.sh quantify  <config.sh> <samples.tsv> <panel_dir> <outdir>
#        all reads mapped to the vOTU panel from the discovery side, counted per vOTU
#   bash run.sh resume <outdir>
#        resubmit what a run did not finish
set -euo pipefail
VP=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
die() { echo "ERROR: $*" >&2; exit 1; }
USAGE="usage: bash run.sh subsample <config.sh> <samples.tsv> <outdir>
       bash run.sh quantify <config.sh> <samples.tsv> <panel_dir> <outdir>
       bash run.sh resume <outdir>"
activate() { set +u; eval "${CONDA_INIT}" || die "CONDA_INIT failed: ${CONDA_INIT}"
             conda activate "${CONDA_ENV}" || die "cannot activate CONDA_ENV=${CONDA_ENV}"; set -u; }

MODE="${1:-}"
case "${MODE}" in
    resume)
        [ $# -eq 2 ] || die "${USAGE}"
        OUT=$(readlink -f "$2")
        [ -s "${OUT}/config.sh" ] && [ -s "${OUT}/manifest.tsv" ] && [ -s "${OUT}/mode" ] || die "${OUT} is not a run directory"
        MODE=$(cat "${OUT}/mode")
        if [ -s "${OUT}/jobs.txt" ]; then
            live=$(squeue -h -o %i -j "$(awk '{print $2}' "${OUT}/jobs.txt" | paste -sd,)" 2> /dev/null || true)
            [ -z "${live}" ] || die "jobs of this run are still queued or running: ${live//$'\n'/ } (wait, or scancel them)"
        fi
        source "${OUT}/config.sh"
        activate
        ;;
    subsample|quantify)
        if [ "${MODE}" = subsample ]; then
            [ $# -eq 4 ] || die "${USAGE}"
            PANEL_SRC=""; OUT=$(readlink -m "$4")
        else
            [ $# -eq 5 ] || die "${USAGE}"
            PANEL_SRC=$(readlink -f "$4"); OUT=$(readlink -m "$5")
        fi
        CONFIG=$(readlink -f "$2"); SAMPLES=$(readlink -f "$3")
        [ -s "${CONFIG}" ] || die "no config file $2"
        [ -s "${SAMPLES}" ] || die "no samplesheet $3"
        [ ! -e "${OUT}" ] || die "${OUT} exists: give a new output directory (or: bash run.sh resume ${OUT})"
        source "${CONFIG}"
        for v in CONDA_INIT CONDA_ENV PARTITION TMPDIR_JOBS MAX_PARALLEL RES_SUBSAMPLE RES_MAP RES_COUNT; do
            [ -n "${!v:-}" ] || die "${v} is empty in ${CONFIG}"
        done
        activate
        if [ "${MODE}" = subsample ]; then tools="seqtk pigz md5sum python3"
        else tools="fastp minimap2 samtools coverm python3"; fi
        for t in ${tools}; do command -v "${t}" > /dev/null || die "${t} not found in ${CONDA_ENV} (environment.yml)"; done
        [ "${MODE}" = subsample ] || python3 -c 'import matplotlib, numpy' 2> /dev/null || die "matplotlib not found in ${CONDA_ENV}"

        mkdir -p "${OUT}/logs" "${OUT}/done"
        python3 -I "${VP}/bin/check_inputs.py" "${SAMPLES}" "${OUT}/manifest.tsv" ${PANEL_SRC:+"${PANEL_SRC}"} || { rm -rf "${OUT}"; exit 1; }
        cp "${SAMPLES}" "${OUT}/samples.tsv"
        cp "${CONFIG}" "${OUT}/config.sh"
        echo "${MODE}" > "${OUT}/mode"
        if [ "${MODE}" = subsample ]; then
            # the samplesheet of the subsampled reads, paths relative to OUT: what the discovery side reads
            awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; print; next}
                {$c["fastq_1"] = "reads/" $c["sample"] "_R1.fq.gz"; $c["fastq_2"] = "reads/" $c["sample"] "_R2.fq.gz"; print}' \
                "${SAMPLES}" > "${OUT}/samples_2M.tsv"
        else
            mkdir -p "${OUT}/panel"
            cp "${PANEL_SRC}/panel.fa" "${PANEL_SRC}/virus_units.tsv" "${PANEL_SRC}/quality_summary.tsv" "${OUT}/panel/"
        fi
        { echo "date: $(date '+%F %T')"; echo "user: ${USER}@$(hostname)"
          echo "command: bash run.sh $*"
          echo "pipeline: ${VP} ($(git -C "${VP}" describe --always --dirty 2> /dev/null || echo 'not a git checkout'))"
          if [ -n "${PANEL_SRC}" ]; then echo "panel: ${PANEL_SRC}"; (cd "${OUT}/panel" && md5sum panel.fa); fi
        } > "${OUT}/run_info.txt"
        conda list > "${OUT}/env_packages.txt"
        ;;
    *) die "${USAGE}" ;;
esac

# --- submission ----------------------------------------------------------------------------------
N=$(wc -l < "${OUT}/manifest.tsv")
sub() {   # name "cpus mem time" "dependency job ids" array-spec script
    local name=$1 dep=$3 arr=$4 script=$5 c m t id
    read -r c m t <<< "$2"
    local opts=(--parsable --kill-on-invalid-dep=yes --partition="${PARTITION}" --cpus-per-task="${c}"
                --mem="${m}" --time="${t}" --job-name="vq_${name}" --chdir="${OUT}"
                --export="ALL,VP=${VP},OUT=${OUT}")
    if [ -n "${arr}" ]; then opts+=(--array="${arr}" --output="${OUT}/logs/${name}-%A_%a.log")
    else opts+=(--output="${OUT}/logs/${name}-%j.log"); fi
    [ -n "${dep}" ] && opts+=(--dependency="afterok:${dep}")
    id=$(sbatch "${opts[@]}" ${SBATCH_EXTRA:-} "${script}") || die "sbatch failed for ${name}"
    id=${id%%;*}
    echo "${name} ${id} $(date '+%F %T')" >> "${OUT}/jobs.txt"
    echo "${id}"
}
all_done() { local s; while IFS=$'\t' read -r s _; do [ -e "${OUT}/done/$1.${s}" ] || return 1; done < "${OUT}/manifest.tsv"; }

S="${VP}/steps"
submitted=""
if [ "${MODE}" = subsample ]; then
    all_done subsample || submitted=$(sub subsample "${RES_SUBSAMPLE}" "" "1-${N}%${MAX_PARALLEL}" "${S}/subsample.sbatch")
    next="when it has finished: send the whole directory ${OUT} (reads/, samples_2M.tsv; check: cd ${OUT}/reads && cat *.md5 | md5sum -c)"
else
    jmap=""
    all_done map || jmap=$(sub map "${RES_MAP}" "" "1-${N}%${MAX_PARALLEL}" "${S}/map.sbatch")
    [ -e "${OUT}/done/count" ] || submitted="${jmap} $(sub count "${RES_COUNT}" "${jmap}" "" "${S}/count.sbatch")"
    next="results: ${OUT}/results/ (viral_rpm.tsv, viral_rpm.pdf, votus.tsv, viral_long.tsv)"
fi
submitted=$(echo ${submitted})
[ -n "${submitted}" ] || { echo "Nothing to submit: every step of ${OUT} is done."; exit 0; }
echo "${MODE}: ${N} samples -> ${OUT}"
echo "submitted jobs: ${submitted}  (list: ${OUT}/jobs.txt; logs: ${OUT}/logs/)"
echo "follow with:    squeue -u ${USER} -n vq_subsample,vq_map,vq_count"
echo "if a job fails: read its log, fix the cause, then  bash ${VP}/run.sh resume ${OUT}"
echo "${next}"
