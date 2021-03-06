#!/usr/bin/env python
import os 
import string
from Bio import SeqIO

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

infile = open('env.fa')
prifile = open('pri.fasta')
outfile = open('ref.fasta','w')
offset  = open('offset','w')
for record in SeqIO.parse(infile,'fasta'):
  refseq = str(record.seq)
prihandle = SeqIO.parse(prifile,'fasta')
for record in prihandle:
  nameF = str(record.id)
  seqF = str(record.seq)
  recordR = prihandle.next()
  nameR = str(recordR.id)
  seqR = str(recordR.seq)
  assert nameF[0:-1] == nameR[0:-1]
  outfile.write('>'+nameF[0:-1]+'\n')
  seqR = rc(seqR)
  index1 = refseq.find(seqF)+len(seqF)
  index2 = refseq.find(seqR)
  amplicon = refseq[index1:index2]
  offset.write(nameF[0:-1]+'\t'+str(index1)+'\t'+str(index2-1)+'\n')
  outfile.write(amplicon+'\n')

infile.close()
prifile.close()
outfile.close()
offset.close()
