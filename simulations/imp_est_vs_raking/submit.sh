#!/bin/bash
# Job name:
#SBATCH --job-name=impute
#
# Partition:
#SBATCH --partition=savio3
#
#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

module load r

R CMD BATCH --no-save run_Y_simple_binomial_Delta_large.R logs/run_Y_simple_binomial_Delta_large.Rout &
R CMD BATCH --no-save run_Y_simple_binomial_Delta_med_strong.R logs/run_Y_simple_binomial_Delta_med_strong.Rout &
R CMD BATCH --no-save run_Y_simple_binomial_Delta_med.R logs/run_Y_simple_binomial_Delta_med.Rout &

wait
