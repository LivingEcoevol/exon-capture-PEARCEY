#!/usr/bin/python3

import os ,sys
from Bio import SeqIO
import glob

path_aa='orthologs/aa/renamed'
path_nt='orthologs/nt/renamed'

if os.path.exists(path_aa):
    pass
else:
    os.makedirs(path_aa)


if os.path.exists(path_nt):
    pass
else:
    os.makedirs(path_nt)


filenames=glob.glob('orthologs/aa/*.fa')
print (len(filenames))

for filename in filenames:
    filename1=filename.replace('orthologs/aa/','')
    gene=filename1.replace('.aa.summarized.fa','')
    out=open('orthologs/aa/renamed/'+filename1,'w')
    for record in SeqIO.parse(filename,'fasta'):
        record.id=record.id.split('|')[1]
        record.description=''
        SeqIO.write(record,out,'fasta')
    out.close()
print ('aa header renamed')

filenames=glob.glob('orthologs/nt/*.fa')
print (len(filenames))

for filename in filenames:
    filename1=filename.replace('orthologs/nt/','')
    gene=filename1.replace('.nt.summarized.fa','')
    out=open('orthologs/nt/renamed/'+filename1,'w')
    for record in SeqIO.parse(filename,'fasta'):
        record.id=record.id.split('|')[1]
        record.description=''
        SeqIO.write(record,out,'fasta')
    out.close()
print ('nt header renamed')
