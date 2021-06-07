#!/usr/bin/python3

import os,sys
import glob
from Bio import SeqIO

# amino acids
path_in='aa_taxa'
path_out='aa_loci'
if os.path.exists(path_out):
    pass
else:
    os.makedirs(path_out)

gene_names=[]
filenames=glob.glob(path_in+'/'+'*.aa.fas')
for filename in filenames:
    for record in SeqIO.parse(filename,'fasta'):
        gene=record.id
        if gene not in gene_names:
            gene_names.append(gene)
print (len(gene_names))

for gene in gene_names:
    outname=gene+'.aa.fas'
    out=open(path_out+'/'+outname,'w')
    for filename in filenames:
        taxon=filename.split('/')[1].split('.')[0]
        for r in SeqIO.parse(filename,'fasta'):
            if gene == r.id:
                r.id=taxon
                r.description=''
                SeqIO.write(r,out,'fasta')
    out.close()
print ('ok')


# nucleaotide
path_in='nt_taxa'
path_out='nt_loci'
if os.path.exists(path_out):
    pass
else:
    os.makedirs(path_out)

gene_names=[]
filenames=glob.glob(path_in+'/'+'*.nt.fas')
for filename in filenames:
    for record in SeqIO.parse(filename,'fasta'):
        gene=record.id
        if gene not in gene_names:
            gene_names.append(gene)
print (len(gene_names))

for gene in gene_names:
    outname=gene+'.nt.fas'
    out=open(path_out+'/'+outname,'w')
    for filename in filenames:
        taxon=filename.split('/')[1].split('.')[0]
        for r in SeqIO.parse(filename,'fasta'):
            if gene == r.id:
                r.id=taxon
                r.description=''
                SeqIO.write(r,out,'fasta')
    out.close()
print ('ok')