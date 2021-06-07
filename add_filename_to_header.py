#!/usr/bin/python3

import os ,sys
from Bio import SeqIO
import glob


filenames=glob.glob('*.fasta')
print (len(filenames))

for filename in filenames:
    taxon=filename.replace('_subseq_coi.fasta','')
    outname=filename.replace('.fasta','_renamed.fas')
    out=open(outname,'w')
    for record in SeqIO.parse(filename,'fasta'):
        record.id=taxon +'_' + record.id
        record.description=''
        SeqIO.write(record,out,'fasta')
    out.close()
    