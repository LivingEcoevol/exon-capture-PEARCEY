#!/usr/bin/python3
#-*- coding: UTF-8 -*-
import os,sys
import glob
from Bio import SeqIO

# amino acids
path_in='orthologs/aa/renamed/'
path_out='aa_taxa'
if os.path.exists(path_out):
    pass
else:
    os.makedirs(path_out)

taxa_names=[]
filenames=glob.glob(path_in+'/'+'*.fa')
for filename in filenames:
    for record in SeqIO.parse(filename,'fasta'):
        taxa=record.id
        if taxa not in taxa_names:
            taxa_names.append(taxa)
print (len(taxa_names))

for taxa in taxa_names:
    outname=taxa+'.aa.fas'
    out=open(path_out+'/'+outname,'w')
    for filename in filenames:
        gene=filename.split('/')[3].split('.')[0]
        for r in SeqIO.parse(filename,'fasta'):
            if taxa == r.id:
                r.id=gene
                r.description=''
                SeqIO.write(r,out,'fasta')
    out.close()
print ('ok')

# nucleotide
path_in='orthologs/nt/renamed/'
path_out='nt_taxa'
if os.path.exists(path_out):
    pass
else:
    os.makedirs(path_out)

taxa_names=[]
filenames=glob.glob(path_in+'/'+'*.fa')
for filename in filenames:
    for record in SeqIO.parse(filename,'fasta'):
        taxa=record.id
        if taxa not in taxa_names:
            taxa_names.append(taxa)
print (len(taxa_names))

for taxa in taxa_names:
    outname=taxa+'.nt.fas'
    out=open(path_out+'/'+outname,'w')
    for filename in filenames:
        gene=filename.split('/')[3].split('.')[0]
        for r in SeqIO.parse(filename,'fasta'):
            if taxa == r.id:
                r.id=gene
                r.description=''
                SeqIO.write(r,out,'fasta')
    out.close()
print ('ok')