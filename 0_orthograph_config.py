#! /usr/bin/env python
import os ,sys
import glob

#####generate config files, given a list of taxa names########
#put this python file into the assembly folder
filenames=glob.glob('*.fasta')

path_config='/scratch1/li266/WGS/WGS_Novogene_2019/orthograph/config' #path to configuration files

#config for relaxed hmmer and blast search, used for orthograph-analyzer
for filename in filenames:
    species=filename.replace('.fasta','')
    print (species)
    input=open(path_config+'/template_analyzer.conf','r')
    output=open(path_config+'/'+species+'_analyzer.conf','w')
    for line in input:
        output.write(line.replace ('template',species))
    input.close()
    output.close()


#config for strict reciprocal search, used for orthograph-reporter
for filename in filenames:
    species=filename.replace('.fasta','')
    print (species)
    input=open(path_config+'/template_reporter.conf','r')
    output=open(path_config+'/'+species+'_reporter.conf','w')
    for line in input:
        output.write(line.replace ('template',species))
    input.close()
    output.close()