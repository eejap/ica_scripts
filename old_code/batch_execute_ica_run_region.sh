#!/bin/bash
#SBATCH -p par-single
#SBATCH -o sbatch_logs/ica_run_regions.out
#SBATCH -e sbatch_logs/ica_run_regions.err
#SBATCH -t 2:00:00  # Adjust the time limit as needed
#SBATCH --mem 60GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6

# Path to your shell script
SHELL_SCRIPT=execute_ica_run_region.sh

# Run your shell script
./$SHELL_SCRIPT sbatch_logs/ica_run_combined.log 2>&1
