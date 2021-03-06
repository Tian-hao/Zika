#!/usr/bin/env python
import sys
from Bio import SeqIO

refdict = {}
reffile = open('../ref/tag')
for line in reffile:
  line = line.rstrip().rsplit('\t')
  refdict[line[0]] = line[1]
reffile.close()
infile1 = open('/u/scratch/y/yushendu/Zika/data/Undetermined_S0_L003_R1_001.fastq')
infile2 = open('/u/scratch/y/yushendu/Zika/data/Undetermined_S0_L003_R2_001.fastq')
outfiles1 = {}
outfiles2 = {}
for bc in refdict:
  outfiles1[bc] = open('/u/scratch/y/yushendu/Zika/split/'+refdict[bc]+'_R1_L2.fastq','w')
  outfiles2[bc] = open('/u/scratch/y/yushendu/Zika/split/'+refdict[bc]+'_R2_L2.fastq','w')
inhandle1 = SeqIO.parse(infile1,'fastq')
inhandle2 = SeqIO.parse(infile2,'fastq')
readcount = 0
for record1 in inhandle1:
  readcount += 1
  if readcount % 100000 == 0: print 'record '+str(readcount)
  record2 = inhandle2.next()
  bc1     = str(record1.seq)[:6]
  bc2     = str(record2.seq)[:6]
  if 'N' in bc1 and 'N' in bc2:
    continue
  elif 'N' in bc1:
    bc = bc2
  elif 'N' in bc2:
    bc = bc1
  elif bc1 != bc2:
    continue
  else:
    bc = bc1
  if bc not in outfiles1: 
    continue
  SeqIO.write(record1,outfiles1[bc],'fastq')
  SeqIO.write(record2,outfiles2[bc],'fastq')
infile1.close()
infile2.close()
for bc in refdict:
  outfiles1[bc].close()
  outfiles2[bc].close()


  
