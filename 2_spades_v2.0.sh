#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=60GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name=spades
#SBATCH --mail-user=living.li@csiro.au

module load spades/3.13.1 jemalloc

export OMP_NUM_THREADS=${SLURM_NTASKS}

assembly_dir="/scratch1/li266/AHE/assembly"
trim_dir="/scratch1/li266/AHE/trim"
threads=10
memory=60

SAMPLES=( $(cat sample_names.txt) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    spades.py \
    --threads ${threads} \
    --memory ${memory} \
    --pe1-1 ${trim_dir}/${SAMPLES[$i]}_R1_paired.fq.gz \
    --pe1-2 ${trim_dir}/${SAMPLES[$i]}_R2_paired.fq.gz \
    --pe1-s ${trim_dir}/${SAMPLES[$i]}_R1_unpaired.fq.gz \
    --pe1-s ${trim_dir}/${SAMPLES[$i]}_R2_unpaired.fq.gz \
    -o ${assembly_dir}/${SAMPLES[$i]}_spades
    cp ${assembly_dir}/${SAMPLES[$i]}_spades/contigs.fasta ${assembly_dir}/${SAMPLES[$i]}_spades.fasta
    rm -r ${assembly_dir}/${SAMPLES[$i]}_spades
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
