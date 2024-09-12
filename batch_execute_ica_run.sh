#!/bin/bash
#SBATCH -p par-single
#SBATCH -o sbatch_logs/ica_run_2_comp.out
#SBATCH -e sbatch_logs/ica_run_2_comp.err
#SBATCH -t 1:00:00  # Adjust the time limit as needed
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6

# Path to your shell script
SHELL_SCRIPT=execute_ica_run.sh

# Run your shell script
./$SHELL_SCRIPT sbatch_logs/ica_run_combined.log 2>&1
