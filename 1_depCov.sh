#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name=dep-cov
#SBATCH --mail-user=living.li@csiro.au

# this scripts is used to calculate the depth of coverage for exon-capture data
# Living 22-Aug-2020

code_dir="/scratch1/li266/program/tccp" #code directory
ref_dir="/scratch1/li266/AHE/orthograph/nt_taxa"  #path to reference
trim_dir="/scratch1/li266/AHE/trim" #path to trimmed reads
#cov_dir="/scratch1/li266/AHE/coverage" ##path to depth of coverage results

module load samtools
module load bwa
module load mosdepth

export OMP_NUM_THREADS=${SLURM_NTASKS}

SAMPLES=( $(cat sample_names.txt) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    
    # copy the reference exon to the working directory
    # extension need to be .fa
    cp ${ref_dir}/${SAMPLES[$i]}_spades.nt.fa .
    ref=${SAMPLES[$i]}_spades.nt.fa
    rgbam=${SAMPLES[$i]}.ReadGrouped.bam
    
    # trimmed reads
    dd1file=${trim_dir}/${SAMPLES[$i]}_R1_trim.fq.gz
    dd2file=${trim_dir}/${SAMPLES[$i]}_R2_trim.fq.gz
    
    # index reference
    bwa index ${ref}
    
    # map reads to the ortholog reference
    bwa mem ${ref} ${dd1file} ${dd2file} | samtools sort -O bam -o ${SAMPLES[$i]}.aln.bam
    samtools index ${SAMPLES[$i]}.aln.bam
    
    #calculate depth of coverage
    mosdepth -n --fast-mode --by 500 ${SAMPLES[$i]}.depcov ${SAMPLES[$i]}.aln.bam

    rm ${SAMPLES[$i]}_spades.nt*

else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
