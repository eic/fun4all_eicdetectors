#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=20:00
#SBATCH --array=0-99
shifter ./AllSi_shifter_hadron.sh $SLURM_ARRAY_TASK_ID 100000 90
