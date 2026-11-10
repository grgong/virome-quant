# Sourced by every step (submitted by run.sh with VP = this directory, OUT = the run directory):
# the run's config, the conda env, a scratch directory removed at exit, and helpers.
set -euo pipefail
: "${VP:?run through run.sh}" "${OUT:?run through run.sh}"
source "${OUT}/config.sh"
set +u; eval "${CONDA_INIT}"; conda activate "${CONDA_ENV}"; set -u
mkdir -p "${TMPDIR_JOBS}"
SCR=$(mktemp -d "${TMPDIR_JOBS}/virome.XXXXXX")
trap 'rm -rf "${SCR}"' EXIT
export TMPDIR="${SCR}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# sample N of manifest.tsv (sample, comma-separated R1 files, comma-separated R2 files)
sample_row() { awk -F'\t' -v n="$1" 'NR==n' "${OUT}/manifest.tsv"; }
log() { echo "[$(date '+%F %T')] $*"; }
