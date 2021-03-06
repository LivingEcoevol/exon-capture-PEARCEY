#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=5Gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=alicut_concat

# Pipeline: aliscore-alicut-concatenation-mare
# Yun (Living) Li version 6-May-2020
# Notes: make sure gene_names.txt in ${base_dir}
# if error: can not find Aliscore module in @INC, then
# specify the path to the module in the commandline: export PERL5LIB=/scratch1/li266/program/alicut


code_dir="/scratch1/li266/program/alicut"
base_dir="/scratch1/li266/AHE/alignment"

mkdir -p ${base_dir}/aa_alicut
mkdir -p ${base_dir}/nt_alicut
mkdir -p ${base_dir}/aa_final
mkdir -p ${base_dir}/nt_final

aa_align_dir=${base_dir}/aa_align
nt_align_dir=${base_dir}/nt_align
aa_alicut_dir=${base_dir}/aa_alicut
nt_alicut_dir=${base_dir}/nt_alicut
aa_final_dir=${base_dir}/aa_final
nt_final_dir=${base_dir}/nt_final

# Fasta to sequential
cd ${aa_align_dir}
for i in *.aa.linsi.fas; do
sample=$(basename $i .aa.linsi.fas)
perl ${code_dir}/fastasingleline.pl $sample.aa.linsi.fas > $sample.aa.ni.linsi.fas
done

cd ${nt_align_dir}
for i in *.nt.linsi.fas; do
sample=$(basename $i .nt.linsi.fas)
perl ${code_dir}/fastasingleline.pl $sample.nt.linsi.fas > $sample.nt.ni.linsi.fas
done

#replace 'forbidden sequence code' generated by MESPA: 'O' -> '-', but be aware that 'O' might be in the headers
#sed -i 's/O/-/g' *.aa.ni.linsi.fas

# Aliscore aa alignment, make sure headers meet Aliscore standard
cd ${aa_align_dir}
for i in *.aa.ni.linsi.fas; do
perl ${code_dir}/Aliscore.02.2.pl -e -i $i
done

# Alicut aa alignment
perl ${code_dir}/ALICUT_V2.3_mod2015.pl -s
mv ALICUT_*.fas ${aa_alicut_dir}

# Copy aa alignment list_of_alicut_sites to nt alignment folder
for i in *.aa.ni.linsi.fas_List_random.txt;do
sample=$(basename $i .aa.ni.linsi.fas_List_random.txt)
cp $sample.aa.ni.linsi.fas_List_random.txt ${nt_align_dir}/$sample.nt.ni.linsi.fas_List_random.txt
done

# Alicut nt alignment based on aa alicut list, does this work when there is intron?
# Double check any error files from the nt alicut results
cd ${nt_align_dir}
perl ${code_dir}/ALICUT_V2.3_mod2015.pl -c -s
mv ALICUT_*.fas ${nt_alicut_dir}

# List genes that have not been cut
cd ${aa_alicut_dir}
ls *.fas >${base_dir}/alicut_gene_names.txt

cd ${base_dir}
sed -i 's/.aa.ni.linsi.fas//g;s/ALICUT_//g' alicut_gene_names.txt

# AA: move Alicut alignment and non-alicut alignment into the final alignment folder
for gene in $(cat gene_names.txt);do
    if grep -qF ${gene} alicut_gene_names.txt;then
        cp ${aa_alicut_dir}/ALICUT_${gene}.aa.ni.linsi.fas ${aa_final_dir}/${gene}.aa.ni.linsi.fas
    else
        cp ${aa_align_dir}/${gene}.aa.ni.linsi.fas ${aa_final_dir}/${gene}.aa.ni.linsi.fas
    fi
done

# NT: move Alicut alignment and non-alicut alignment into the final alignment folder
for gene in $(cat gene_names.txt);do
    if grep -qF ${gene} alicut_gene_names.txt;then
        cp ${nt_alicut_dir}/ALICUT_codon_${gene}.nt.ni.linsi.fas ${nt_final_dir}/${gene}.nt.ni.linsi.fas
    else
        cp ${nt_align_dir}/${gene}.nt.ni.linsi.fas ${nt_final_dir}/${gene}.nt.ni.linsi.fas
    fi
done

# Concatenation
cd ${nt_final_dir}
perl ${code_dir}/FASconCAT-G_v1.04.pl -s -l

cd ${aa_final_dir}
perl ${code_dir}/FASconCAT-G_v1.04.pl -s -l

# MARE subsetting the aa matrix
module load mare

# Modify aa partition file
sed -i 's/LG, /charset /g;s/.aa.ni.linsi//g;s/$/;/g' FcC_supermatrix_partition.txt

# Run mare on aa matrix
MARE FcC_supermatrix_partition.txt FcC_supermatrix.fas -d 2 -t 1.5 -m


