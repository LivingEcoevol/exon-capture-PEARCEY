#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=1GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=orthograph

result_dir="/datasets/work/ncmi-cwgs/work/temp_orthograph_LY/orthograph_results"

mkdir -p orthologs
perl /apps/orthograph/0.6.3/summarize_orthograph_results.pl -i ${result_dir} -o orthologs -c -t -s 
#-d species_to_exclude.txt

mv orthologs/aa_summarized orthologs/aa
mv orthologs/nt_summarized orthologs/nt


