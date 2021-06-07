#!/usr/bin/python3

import os ,sys
from Bio import SeqIO
import glob

name_list='name_list.txt'
dic={}
for i in open(name_list).readlines()[1:]:
    key=i.split('\t')[0].strip()
    value=i.split('\t')[1].strip()
    dic[key]=value
print (dic)

filenames=glob.glob('*.fas_reduced')
print (len(filenames))

for filename in filenames:
    outname=filename.replace('.fas_reduced','.renamed.fa')
    out=open(outname,'w')
    for record in SeqIO.parse(filename,'fasta'):
        #modify the pattern case by case
        voucher=record.id.replace('_spades','')
        if voucher in dic:
            record.id=dic[voucher]
            record.description=''
            SeqIO.write(record,out,'fasta')
        else:
            SeqIO.write(record,out,'fasta')
    out.close()
    