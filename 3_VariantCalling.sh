#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name=TCCP_VC
#SBATCH --mail-user=living.li@csiro.au

# this script aims to call SNPs from exon-capture data
# Living 22-Aug-2020

code_dir="/scratch1/li266/program/tccp" #code directory
ref_dir="/scratch1/li266/AHE/orthograph/nt_taxa"  #path to exon reference
trim_dir="/scratch1/li266/AHE/trim" #path to trimmed reads

module load samtools
module load bwa
module load bcftools
module load vcftools

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
    
    # map reads to the exon reference
    bwa mem ${ref} ${dd1file} ${dd2file} | samtools sort -O bam -o ${SAMPLES[$i]}.aln.bam
    samtools index ${SAMPLES[$i]}.aln.bam

    #remove unmapped reads from bam file
    samtools view -b -F 4 ${SAMPLES[$i]}.aln.bam > ${SAMPLES[$i]}.mapped.bam

    # prepare gatk, read grouping
    java -jar ${code_dir}/AddOrReplaceReadGroups.jar INPUT=${SAMPLES[$i]}.mapped.bam OUTPUT=${rgbam} RGID=${SAMPLES[$i]} RGLB=${SAMPLES[$i]} RGPU=${SAMPLES[$i]} RGPL=ILLUMINA RGSM=${SAMPLES[$i]}
    # index reference
    samtools faidx ${ref}
    # index reading group
    samtools index ${rgbam}
    # make dictionary
    java -jar ${code_dir}/CreateSequenceDictionary.jar R=${ref} O=${ref/fa/dict}

    # call haplotypes
    java -jar ${code_dir}/GenomeAnalysisTK.jar -R ${ref} -I ${rgbam} -T HaplotypeCaller -o ${SAMPLES[$i]}.vcf
    # perform read backed phasing
    java -jar ${code_dir}/GenomeAnalysisTK.jar -R ${ref} -I ${rgbam} -T ReadBackedPhasing --variant ${SAMPLES[$i]}.vcf --min_base_quality_score 21 -o ${SAMPLES[$i]}.ReadBackedPhased.vcf
    
    # remove indels outright
    vcftools --vcf ${SAMPLES[$i]}.ReadBackedPhased.vcf --remove-indels --recode --recode-INFO-all --out ${SAMPLES[$i]}.ReadBackedPhased
    
    # calculate coverage
    java -jar ${code_dir}/GenomeAnalysisTK.jar -R ${ref} -I ${rgbam} -T DepthOfCoverage --omitIntervalStatistics -o ${SAMPLES[$i]}.DepthOfCoverageTable

    # process vcf file
    bgzip -c ${SAMPLES[$i]}.ReadBackedPhased.recode.vcf > ${SAMPLES[$i]}.vcf.gz
    bcftools index ${SAMPLES[$i]}.vcf.gz

    # haplotypes
    bcftools consensus -H 1 -f ${ref} -o ${SAMPLES[$i]}_H1.fas ${SAMPLES[$i]}.vcf.gz
    bcftools consensus -H 2 -f ${ref} -o ${SAMPLES[$i]}_H2.fas ${SAMPLES[$i]}.vcf.gz
    
    # diplotypes
    bcftools consensus -i -f ${ref} -o ${SAMPLES[$i]}_D.fas ${SAMPLES[$i]}.vcf.gz

    # mask out sequences below a minimum coverage threshold; f: the fasta file; c: the coverage file; o: the coverage cutoff
    python ${code_dir}/mask_low_cov2.py -f ${SAMPLES[$i]}_H1.fas -c ${SAMPLES[$i]}.DepthOfCoverageTable -o 10 > ${SAMPLES[$i]}_masked_H1.fas
    python ${code_dir}/mask_low_cov2.py -f ${SAMPLES[$i]}_H2.fas -c ${SAMPLES[$i]}.DepthOfCoverageTable -o 10 > ${SAMPLES[$i]}_masked_H2.fas
    python ${code_dir}/mask_low_cov2.py -f ${SAMPLES[$i]}_D.fas -c ${SAMPLES[$i]}.DepthOfCoverageTable -o 10 > ${SAMPLES[$i]}_masked_D.fas
    
    # remove unnecessary files
    ${SAMPLES[$i]}.DepthOfCoverageTable*
    ${SAMPLES[$i]}.ReadBackedPhased*
    ${SAMPLES[$i]}.aln.bam*
    ${SAMPLES[$i]}.mapped.bam*
    ${SAMPLES[$i]}.vcf*
    ${SAMPLES[$i]}_spades.nt*
    
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
