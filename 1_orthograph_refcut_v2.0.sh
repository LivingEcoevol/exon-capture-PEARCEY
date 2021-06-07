#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=orth_ref_cut

#this script is used to make 100% corresponding nt and aa reference sequences
# Yun (Living) Li 6-May-2020

#remove gaps
sed -i 's/-//g' aa_clean/*.fas
sed -i 's/-//g' nt_clean/*.fas

mkdir -p ref_cut

for sample in $(cat sample_names.txt);do
perl /apps/orthograph/0.6.3/make-ogs-corresponding.pl --num-threads 10 --outdir ref_cut aa_clean/$sample.ni.aa.fas nt_clean/$sample.ni.nt.fas
done
