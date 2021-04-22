#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=1:00:00
#SBATCH --array=0-199
shifter ./AllSi_shifter.sh $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID 10000 0.3 0.3 0.3 0. 30. 3
