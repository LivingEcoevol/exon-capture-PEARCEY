#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=500Mb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=orthograph

module load orthograph

#modify config files to use strict search on Tribolium_castaneum
config_dir="/scratch1/li266/AHE/orthograph/config"

SAMPLES=( $(cat sample_names.txt) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
	i=$SLURM_ARRAY_TASK_ID
	orthograph-reporter -c ${config_dir}/${SAMPLES[$i]}_reporter.conf
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi