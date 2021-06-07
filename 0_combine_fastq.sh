#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=merge_fastq

for sample in $(cat sample_names.txt); do
cat ${sample}_R1_paired.fq.gz ${sample}_R1_unpaired.fq.gz > ${sample}_R1_trim.fq.gz
cat ${sample}_R2_paired.fq.gz ${sample}_R2_unpaired.fq.gz > ${sample}_R2_trim.fq.gz
echo ${sample} "done"
done