#!/bin/bash
## Project:
#SBATCH --account=nn9998k
## Job name:
#SBATCH --job-name=myjob
## Wall time limit:
#SBATCH --time=1-0:0:0
## Number of nodes:
#SBATCH --nodes=4
## Number of tasks to start on each node:
#SBATCH --ntasks-per-node=1
## Set OMP_NUM_THREADS
#SBATCH --cpus-per-task=1

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default
module load R/3.4.0-intel-2017a-X11-20170314
module list             # List loaded modules, for easier debugging

## Prepare input files
cp -r ../../EFSIS $SCRATCH
cd $SCRATCH/EFSIS/scripts

## Make sure output is copied back after job finishes
savefile auc.all*

## Run the application
srun ./efsis_diff_num_sel.R

#export TMP=/cluster/work/users/xzh004/
#export TMPDIR=$TMP
#
#aprun -n1 -N1 -m30000M R --no-save < ./efsis_diff_num_sel.R
