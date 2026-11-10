# Virome quantification: copy this file, fill it in, and pass it to run.sh.
# A copy is kept in every output directory.

# --- software -------------------------------------------------------------------------------
CONDA_INIT=". /local/env/envconda.sh"      # how a job makes `conda` available (GenOuest); elsewhere e.g. "source ~/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV=/path/to/envs/virome_quant       # the env created from environment.yml (path or name)

# --- cluster --------------------------------------------------------------------------------
PARTITION=genouest                         # SLURM partition (`sinfo -s` lists them)
SBATCH_EXTRA=""                            # any other sbatch options, e.g. "--account=myproject"
TMPDIR_JOBS=/scratch/$USER/tmp             # temporary files of the jobs: up to ~20 GB per running mapping job
MAX_PARALLEL=20                            # array tasks running at the same time (bounds TMPDIR_JOBS use)

# resources per step: cpus memory time
RES_SUBSAMPLE="2 8G 02:00:00"
RES_MAP="8 16G 04:00:00"
RES_COUNT="4 16G 01:00:00"
