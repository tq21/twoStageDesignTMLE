#!/bin/bash
# Job name:
#SBATCH --job-name=ipcw_pi
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
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=sky.qiu@berkeley.edu

module load r

# R CMD BATCH --no-save run_ver_0.R logs/run_ver_0.Rout &
# R CMD BATCH --no-save run_ver_1.R logs/run_ver_1.Rout &
# R CMD BATCH --no-save run_ver_2.R logs/run_ver_2.Rout &
# R CMD BATCH --no-save run_ver_3.R logs/run_ver_3.Rout &
R CMD BATCH --no-save run_plugin_ver_0.R logs/run_plugin_ver_0.Rout &
wait
