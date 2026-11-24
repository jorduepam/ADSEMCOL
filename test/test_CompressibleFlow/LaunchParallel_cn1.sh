#!/bin/bash

# https://guiesbibtic.upf.edu/recerca/hpc/array-jobs

# NOTE: SBATCH parameters must be defined first. Do not define any variable before them.

#SBATCH --job-name=LP3
#SBATCH --partition=cn1
#SBATCH --array=1-10%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=30-00:00:00
#SBATCH --output=temp/%x_%A_%a.out
#SBATCH --error=temp/%x_%A_%a.err

# Launch independent simulations as a Slurm job array. From the repository root:
#
#   mkdir -p temp
#   sbatch test/test_CompressibleFlow/LaunchParallel_cn1.sh
#
# Slurm opens the output files before the script runs, so temp must already exist.
# The job name selects test/test_CompressibleFlow/<job-name>.lp: task 1 reads
# line 1, and so on. Create the matching file and check the array range before
# submitting. Arguments must be separated by spaces or tabs; quotes, shell
# expansion and spaces within arguments are not supported.
#
# Use SC identifiers distinct from all other simulations, including pending
# ones. A different job name does not prevent mesh or result files from being
# overwritten when the same SC is used.
#
# Candidate configuration from a cn1 benchmark: up to 5 simulations at once,
# with 4 Julia threads and 1 BLAS thread each. This is not necessarily optimal
# for other cases or long simulations.
#
# Change the directives above or pass options to sbatch:
#   --partition=cn1       partition; cn2 can also be selected
#   --array=1-10%5        lines 1 to 10; at most 5 tasks at once
#   --cpus-per-task=4     CPUs per simulation; also Julia threads by default
#   --mem=40G             memory requested per simulation, not for the array
#   --time=30-00:00:00    30 days per task; adjust for the case and policy
#
# Example: sbatch --array=1-4%2 --cpus-per-task=8 --time=02:00:00 \
#              test/test_CompressibleFlow/LaunchParallel_cn1.sh
#
# To set Julia and BLAS threads separately from the allocated CPUs:
#   sbatch --cpus-per-task=8 \
#     --export=ALL,ADSEMCOL_JULIA_THREADS=4,ADSEMCOL_BLAS_THREADS=1 \
#     test/test_CompressibleFlow/LaunchParallel_cn1.sh
# Julia defaults to the allocated CPU count and BLAS to 1 thread. Both must
# be between 1 and the allocated CPU count. These checks do not rule out
# oversubscription: extra BLAS threads can hurt performance; benchmark first.
#
# Before a large array, precompile packages once on a suitable compute node
# using the same Julia module, active environment and depot as the jobs:
#   julia -e 'using Pkg; Pkg.precompile()'
# This reduces concurrent compilation, but does not guarantee that all code is
# compiled or that first-start overhead is eliminated.
#
# Slurm saves this script at submission; pending tasks read the .lp and Julia
# code when they start. Do not modify those files while tasks need them. Log
# names contain the array ID. Every 500 accepted time steps, .out reports
# Julia's peak RSS and, when available, Slurm cgroup memory and limits in GiB.
# Set ADSEMCOL_MEMORY_REPORT_STEPS via --export to change the interval, or
# set it to 0 to disable periodic reports. The small reporting cost is not zero.
# Julia's final peak RSS is also printed in bytes on normal exit, including
# caught errors. After srun ends, the batch shell reports the job cgroup peak
# if available, including on failure. A forced exit can prevent either report;
# the last sample may precede an OOM kill. See the diagnostic for limitations.
#
# cn1: srun has failed to communicate with the controller. Resolve this with
# the administrator and test a small array before using this for production.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
cpus="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is missing}"
julia_threads="${ADSEMCOL_JULIA_THREADS:-$cpus}"
blas_threads="${ADSEMCOL_BLAS_THREADS:-1}"

for value in "$cpus" "$julia_threads" "$blas_threads"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: CPU and thread counts must be positive integers: $value" >&2
        exit 1
    fi
done
if (( julia_threads > cpus || blas_threads > cpus )); then
    echo "Error: Julia or BLAS threads exceed the $cpus allocated CPUs." >&2
    exit 1
fi

input="test/test_CompressibleFlow/${SLURM_JOB_NAME:?SLURM_JOB_NAME is missing}.lp"
if [[ ! -r "$input" ]]; then
    echo "Error: cannot read the arguments file $input" >&2
    exit 1
fi
line=$(sed -n "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is missing}p" "$input")
line=${line%$'\r'}
read -r -a args <<< "$line"
if (( ${#args[@]} == 0 )) || [[ "${args[0]}" == \#* ]]; then
    echo "Error: line ${SLURM_ARRAY_TASK_ID} is empty or commented out in $input" >&2
    exit 1
fi

module load julia
module load bamg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS="$blas_threads" ADSEMCOL_BLAS_THREADS="$blas_threads"
export ADSEMCOL_REPORT_MAXRSS=1
export ADSEMCOL_MEMORY_REPORT_STEPS="${ADSEMCOL_MEMORY_REPORT_STEPS:-500}"

echo "Start: $(date); node: $(hostname); array: $SLURM_ARRAY_JOB_ID; task: $SLURM_ARRAY_TASK_ID"
echo "CPUs: $cpus; Julia threads: $julia_threads; BLAS threads: $blas_threads; input: $input"

report_job_peak() {
    local hierarchy controllers path job_path peak_file peak
    if [[ ! "${SLURM_JOB_ID:-}" =~ ^[0-9]+$ || ! -r /proc/self/cgroup ]]; then
        echo "Job cgroup peak unavailable: missing job ID or cgroup information."
        return
    fi
    while IFS=: read -r hierarchy controllers path; do
        if [[ "$hierarchy" == 0 && -z "$controllers" ]]; then
            path="${path%/}/"
            case "$path" in
                *"/job_${SLURM_JOB_ID}/"*)
                    job_path="${path%%/job_"${SLURM_JOB_ID}"/*}/job_${SLURM_JOB_ID}"
                    peak_file="/sys/fs/cgroup${job_path}/memory.peak"
                    if [[ -r "$peak_file" ]] && peak=$(<"$peak_file") && [[ "$peak" =~ ^[0-9]+$ ]]; then
                        echo "Job cgroup peak after srun: $peak bytes"
                    else
                        echo "Job cgroup peak unavailable: cannot read $peak_file."
                    fi
                    return
                    ;;
            esac
        fi
    done < /proc/self/cgroup
    echo "Job cgroup peak unavailable: no accessible Slurm cgroup v2 found."
}

# One Julia process per task, bound to the CPUs assigned to the Slurm step.
# Additional internal library threads are not necessarily active simultaneously.
status=0
srun --ntasks=1 --cpus-per-task="$cpus" --cpu-bind=cores \
    julia -t "$julia_threads" test/test_CompressibleFlow/Launch.jl "${args[@]}" || status=$?

report_job_peak
if (( status != 0 )); then
    echo "srun failed with exit status $status; inspect Slurm state and the error log." >&2
    exit "$status"
fi
echo "Julia finished: $(date)."