#!/bin/bash
# this script is to automate OGS building process in Orthograph
# Yun (Living) Li 6-May-2020

module load orthograph

# create tables and SQLite database
orthograph-manager --create

# load aa sequences into the database
for sample in $(cat sample_names.txt);do
orthograph-manager --load-ogs-peptide ref_cut/corresp-${sample}.ni.aa.fas --ogs-version 10 --ogs-taxon-name ${sample}
done

# load nt sequences into the database
for sample in $(cat sample_names.txt);do
orthograph-manager --load-ogs-nucleotide ref_cut/corresp-${sample}.ni.nt.fas --ogs-version 10 --ogs-taxon-name ${sample}
done

#upload your ortholog gene set into the database
#make a note of your ortholog set name
orthograph-manager COG_tab.txt

#print a list of all OGS present in the database
orthograph-manager -lo

#print a list of taxa present in the database
orthograph-manager -lt

#check uploaded sets
orthograph-manager -ls

