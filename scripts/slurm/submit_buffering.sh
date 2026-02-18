#!/bin/bash
#SBATCH --job-name=buffering_optim
#SBATCH --output=logs/buffering_%A_%a.out  # %A=JobId, %a=ArrayIdx
#SBATCH --error=logs/buffering_%A_%a.err

#SBATCH --array=1-100       # Example: Run 100 optimization starts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Optim is typically single-threaded unless parallelized internally
#SBATCH --mem=4G            # Adjust based on 'run_sim_jobs_in_memory' needs
#SBATCH --time=00:30:00     
#SBATCH --qos=normal

# ---- Container Config ----
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

# ---- Execution Config ----
ARRAY_IDX=$SLURM_ARRAY_TASK_ID

echo "Starting Buffering Optimization Array Job on $(hostname)"
echo "Array Index: ${ARRAY_IDX}"

# ---- Create Logs Directory ----
mkdir -p logs

# ---- FORCE CACHE REFRESH ----
id $USER > /dev/null 2>&1
sleep 2

# ---- Execution ----
# R Script Args: 
# 1: ARRAY_IDX (Used for seed and output filename)

apptainer exec $BINDS $CONTAINER_URI \
  Rscript scripts/fit_wrappers/fit_buffering.R \
  "$ARRAY_IDX"