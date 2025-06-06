#!/bin/bash
# Job name:
#SBATCH --job-name=twoStage
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
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

module load r

R CMD BATCH --no-save run_Y_1.0_Delta_1.0.R logs/run_Y_1.0_Delta_1.0.Rout &
R CMD BATCH --no-save run_Y_1.0_Delta_1.1.R logs/run_Y_1.0_Delta_1.1.Rout &
R CMD BATCH --no-save run_Y_1.1_Delta_1.0.R logs/run_Y_1.1_Delta_1.0.Rout &
R CMD BATCH --no-save run_Y_1.1_Delta_1.1.R logs/run_Y_1.1_Delta_1.1.Rout &

wait
