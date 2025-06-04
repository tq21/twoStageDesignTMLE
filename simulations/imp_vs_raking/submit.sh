#!/bin/bash
# Job name:
#SBATCH --job-name=stupid
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
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=10
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

module load r

R CMD BATCH --no-save run_miss_type_med_1000.R logs/run_miss_type_med_1000.Rout &
R CMD BATCH --no-save run_miss_type_med_2000.R logs/run_miss_type_med_2000.Rout &
R CMD BATCH --no-save run_miss_type_med_3000.R logs/run_miss_type_med_3000.Rout &
R CMD BATCH --no-save run_miss_type_med_4000.R logs/run_miss_type_med_4000.Rout &

wait
