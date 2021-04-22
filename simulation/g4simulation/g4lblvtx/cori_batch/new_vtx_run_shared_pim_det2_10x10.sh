#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=1:00:00
#SBATCH --array=0-499
shifter ./AllSi_shifter_specify_det_and_pix.sh $SLURM_ARRAY_TASK_ID 10000 "pi-" 2 10
