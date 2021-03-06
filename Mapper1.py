#!/usr/bin/env python
import os
import sys
import string
import operator
from itertools import imap
from Bio import SeqIO

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def hamming(str1,str2):
  assert len(str1) == len(str2)
  return sum(imap(operator.ne, str1, str2))

def readrefs(infile):
  refdict  = {}
  inhandle = open(infile)
  for record in SeqIO.parse(inhandle,'fasta'):
    name = str(record.id)
    seq  = str(record.seq)
    refdict[name] = seq
  inhandle.close()
  return refdict

def anchor_offset(seq,primers):
  for primer in primers:
    pri_seq = primers[primer]
    for pri_offset in range(0,5):
      pri = pri_seq[pri_offset:]
      for seq_offset in range(0,5):
        pseq = seq[seq_offset:seq_offset+len(pri)]
        if hamming(pseq,pri) < 3:
          direction = primer[-1]
          amp = primer[:-1]
          offset = seq_offset+len(pri)
          return direction, amp, offset
  return 'NA','NA',0

def main():
  workpath = '/u/flashscratch/t/tianhao/Zika/'
  infile1  = workpath+'split/'+sys.argv[1]
  infile2  = infile1.replace('_R1_','_R2_')
  filenum  = sys.argv[1].rsplit('fastq')[0].rsplit('_') 
  filenum  = filenum[2] +'_'+ filenum[-1]
  outfile  = workpath+'Zika/'+filenum+'m'
  refseqs  = readrefs('../ref/ref.fasta')
  primers  = readrefs('../ref/pri.fasta')
  bclength = 6
  bcerror  = 2

  count    = [0,0,0,0,0,0]
  #count[0] adaptors are not mathced.
  #count[1] unmapped reads
  #count[2] aborted reads
  #count[3] strange mapped reads
  #count[4] wrong read pairs
  #count[5] low quality reads
  handle1  = open(infile1)
  handle2  = open(infile2)
  handleo  = open(outfile,'w')
  seqhand1 = SeqIO.parse(handle1,'fastq')
  seqhand2 = SeqIO.parse(handle2,'fastq')
  for readcount,record1 in enumerate(seqhand1):
    if readcount % 10000 == 0:
       print 'processing line '+str(readcount)
    record2 = seqhand2.next()
    assert  record1.id == record2.id
    seq1    = str(record1.seq)
    seq2    = str(record2.seq)
    qual1   = record1.letter_annotations["phred_quality"]
    qual2   = record2.letter_annotations["phred_quality"]
    bc1     = seq1[:bclength]
    bc2     = seq2[:bclength]
    if hamming(bc1,bc2) > bcerror or ('N' in bc1 and 'N' in bc2): 
      count[0] += 1
      #print bc1,bc2
      #print seq1
      #print seq2
      continue
    if 'N' in bc1: bc = bc2
    else: bc = bc1

    #aborted reads checking    
    if len(seq1) < 93 or len(seq2) < 93: 
      count[2] += 1
      continue

    #offset anchoring
    seq1 = seq1[bclength:]
    seq2 = seq2[bclength:]
    direction1, amp1, offset1 = anchor_offset(seq1,primers)
    direction2, amp2, offset2 = anchor_offset(seq2,primers)
    if direction1 == 'NA' or direction2 == 'NA':
      count[1] += 1
      continue
    if amp1 != amp2 or direction1 == direction2:
      count[3] += 1
      continue
    refseq = refseqs[amp1]
    seq1   = seq1[offset1:]
    seq2   = seq2[offset2:]
    qual1  = qual1[bclength+offset1:]
    qual2  = qual2[bclength+offset2:]
    assert len(seq1) == len(qual1)
    assert len(seq2) == len(qual2)
    
    #adjust length
    lseq1 = len(seq1) ; lseq2 = len(seq2) ; lref = len(refseq)
    dlen1 = 0 ; dlen2 = 0
    if lseq1 < lref:
      dlen1  = lref-lseq1
      seq1  += refseq[lseq1:]
      qual1 += dlen1*[100]
    if lseq2 < lref:
      dlen2  = lref-lseq2
      seq2  += refseq[lseq2:]
      qual2 += dlen2*[100]
    if lseq1 > lref:
       seq1  = seq1[:lref]
       qual1 = qual1[:lref]
    if lseq2 > lref:
       seq2  = seq2[:lref]
       qual2 = qual2[:lref]
    if direction1 == 'R': 
      seq1  = rc(seq1)
      qual1 = qual1[::-1]
      begin_pos = dlen1
      end_pos = lref-dlen2
    else:
      seq2  = rc(seq2)
      qual2 = qual2[::-1]
      begin_pos = dlen2
      end_pos = lref-dlen1
    assert len(seq1)  == len(refseq)
    assert len(seq2)  == len(refseq)
    assert len(qual1) == len(refseq)
    assert len(qual2) == len(refseq)

    #call mutation
    muts = []
    low_qs_flag = 0
    unpair_flag = 0
    for n in range(begin_pos,end_pos):
      if seq1[n] != refseq[n] and seq1[n] == seq2[n] and qual1[n] >= 30 and qual2[n] >= 30: 
        mut = refseq[n] + str(n+1) + seq1[n]
        muts.append(mut)
      elif seq1[n] != refseq[n] and seq1[n] == seq2[n] and (qual1[n] < 30 or qual2[n] < 30): 
        low_qs_flag = 1
        #break
      elif seq1[n] != seq2[n]:
        unpair_flag = 1
        #break
    if len(muts) == 0: muts = 'WT'
    else: muts = '-'.join(muts)
    if low_qs_flag == 1: count[5] += 1; 
    if unpair_flag == 1: count[4] += 1;
    handleo.write(amp1+'\t'+bc+'\t'+muts+'\n')

  handleo.write('Adaptors are not mathced: '+str(count[0])+'\n')
  handleo.write('unmapped reads: '+str(count[1])+'\n')
  handleo.write('aborted reads: '+str(count[2])+'\n')
  handleo.write('strange mapped reads: '+str(count[3])+'\n')
  handleo.write('wrong read pairs: '+str(count[4])+'\n')
  handleo.write('low quality reads: '+str(count[5])+'\n')
  
  handle1.close()
  handle2.close()
  handleo.close()

if __name__ == '__main__':
  main()
