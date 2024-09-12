#!/bin/bash
#SBATCH -p par-single
#SBATCH -o sbatch_logs/ica_run_t_4_comp_028A.out
#SBATCH -e sbatch_logs/ica_run_t_4_comp_028A.err
#SBATCH -t 14:00:00  # Adjust the time limit as needed
#SBATCH --mem 50GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6

# Path to your shell script
SHELL_SCRIPT=execute_ica_run_frames.sh

# Run your shell script
./$SHELL_SCRIPT sbatch_logs/ica_run_combined.log 2>&1
