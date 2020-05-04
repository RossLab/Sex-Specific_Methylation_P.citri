#!/usr/bin/env hmarshall

# This creates a new fasta file containing just the sequence of the annotated coordinates of a given gff
# Need to rename the .fa to .fas
# Can then just run on the command line python <script.py>
# Not sure if the below libraries are needed, script taken from Sam Lewis Github

import os
import sys
import subprocess
import random as random
import shutil
import re
import argparse

def Extractor(fastafile,gff):
	# extract sequences corresponding to gff annotations to a fasta file (NB: this forces strandedness i.e. reverse-complements annotations on the antisense strand)
	cmd = 'bedtools getfasta -s -fo ' + fastafile.replace('.fas','_TE.fas') + ' -fi ' + fastafile + ' -bed ' + gff
	subprocess.call(cmd,shell=True)
	print('Fasta file written: ' + fastafile.replace('.fas','_TE.fas'))
	# remove temp/intermediate files
	os.remove(fastafile + '.fai')

Extractor('PCITRI.assembly.v0.fas', 'PCITRI.assembly.v0.fa.out.gff')

# Remove duplicate rows
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' PCITRI.assembly.v0_TE.fasta | \
	awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > PCITRI.assembly.v0_TE_nodups.fasta

# When this is done split the huge fasta file into many smaller ones for interproscan 
pyfasta split -n 40 PCITRI.assembly.v0_TE_nodups.fasta 

