#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=AHE_trim
#SBATCH --mail-user=living.li@csiro.au

# Pipeline: fastqc-deduplication-trim
# Yun (Living) Li, 6-May-2020

#adapter_dir='/datasets/work/ncmi-cwgs/work/' #directory of the bbduk adaptor file
adaptor_file='/datasets/work/ncmi-cwgs/work/TruSeq3-PE-2.fa' #location of the adaptor file for trimmomatic
read_dir="/scratch1/li266/AHE/reads"
trim_dir="/scratch1/li266/AHE/trim"
fastqc_dir="/scratch1/li266/AHE/fastqc"

module load fastqc/0.11.8
module load fastuniq/1.1
#module load bbmap/38.37
module load trimmomatic/0.38


SAMPLES=( $(cat sample_names.txt) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    #remove PCR duplicates
    gunzip -k ${read_dir}/${SAMPLES[$i]}_R1_001.fastq.gz ${read_dir}/${SAMPLES[$i]}_R2_001.fastq.gz
    mv ${read_dir}/${SAMPLES[$i]}_R1_001.fastq .
    mv ${read_dir}/${SAMPLES[$i]}_R2_001.fastq .
    ls ${SAMPLES[$i]}_R1_001.fastq ${SAMPLES[$i]}_R2_001.fastq > ${SAMPLES[$i]}_read_list.txt
    fastuniq -i ${SAMPLES[$i]}_read_list.txt -t q -o ${SAMPLES[$i]}_R1_dedup.fq -p ${SAMPLES[$i]}_R2_dedup.fq

    #pre trimming fastqc
    fastqc ${SAMPLES[$i]}_R1_dedup.fq ${SAMPLES[$i]}_R2_dedup.fq ${SAMPLES[$i]}_R1_001.fastq ${SAMPLES[$i]}_R2_001.fastq  -o ${fastqc_dir}
    rm ${SAMPLES[$i]}_R1_001.fastq ${SAMPLES[$i]}_R2_001.fastq

    #trim reads with bbduk
    #bbduk.sh t=10 in1=${SAMPLES[$i]}_R1_dedup.fq in2=${SAMPLES[$i]}_R2_dedup.fq out1=${trim_dir}/${SAMPLES[$i]}_R1_trim.fq.gz out2=${trim_dir}/${SAMPLES[$i]}_R2_trim.fq.gz overwrite=true qtrim=rl trimq=20 minlen=35 k=23 ktrim=r mink=11 hdist=1 ref=${adapter_dir}/adapters.fa 
    #rm ${SAMPLES[$i]}_R1_dedup.fq ${SAMPLES[$i]}_R2_dedup.fq

    #trim reads with trimmomatic
    trimmomatic PE -phred33 ${SAMPLES[$i]}_R1_dedup.fq ${SAMPLES[$i]}_R2_dedup.fq ${trim_dir}/${SAMPLES[$i]}_R1_paired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R1_unpaired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R2_paired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R2_unpaired.fq.gz ILLUMINACLIP:${adaptor_file}:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:40
    rm ${SAMPLES[$i]}_R1_dedup.fq ${SAMPLES[$i]}_R2_dedup.fq

    #post trimming fastqc
    #fastqc ${trim_dir}/${SAMPLES[$i]}_R1_trim.fq.gz ${trim_dir}/${SAMPLES[$i]}_R2_trim.fq.gz  -o ${fastqc_dir}
    fastqc ${trim_dir}/${SAMPLES[$i]}_R1_paired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R1_unpaired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R2_paired.fq.gz ${trim_dir}/${SAMPLES[$i]}_R2_unpaired.fq.gz  -o ${fastqc_dir}
    
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi