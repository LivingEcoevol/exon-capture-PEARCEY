#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=1GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=orthograph_building

module load orthograph

config_dir="/scratch1/li266/AHE/orthograph/config"

# create the required database structure for the Orthograph search program
orthograph-analyzer -c ${config_dir}/building_sqlite.conf