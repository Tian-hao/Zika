#!/usr/bin/env python
import os
import sys
import glob
import string

dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", "XXX":"X"}

def anchor_mut(mut,amp,ampdict):
  pos = int(mut[1:-1])
  pos += ampdict[amp][0]
  mut = mut[0] + str(pos) + mut[-1]
  return mut

def translate_mut(aa_geno,seq):
  pos      = int(aa_geno[1:-1])
  wtcodon  = seq[pos/3*3:pos/3*3+3]
  seq      = seq[:pos]+aa_geno[-1]+seq[pos+1:]
  mutcodon = seq[pos/3*3:pos/3*3+3]
  if 'N' in mutcodon: mutaa = 'X'
  else: mutaa = dnamap[mutcodon]
  wtaa = dnamap[wtcodon]
  return wtaa+str(pos/3+1)+mutaa
  
def main():
  workpath = '/u/flashscratch/t/tianhao/Zika/Zika/'
  outfile  = open('../result/readcount.txt','w')
  depfile  = open('../result/depth.txt','w')
  sumfile  = open('../result/summary.txt','w')

  #read reference files  
  genodict = {} #genotype: lib_name: count
  libdict  = {} #barcode: lib_name
  depdict  = {} #lib_name: amp: count
  wtdict   = {} #lib_name: amp: count
  sumdict  = {} #summary: count
  liblist  = []
  ORFlist  = []
  tagfile  = open('../ref/tag')
  for line in tagfile:
    line = line.rstrip().rsplit('\t')
    libdict[line[0]] = line[1]
    depdict[line[1]] = {}
    wtdict[line[1]]  = {}
    liblist.append(line[1])
  tagfile.close()

  ORFdict  = {}
  ORFfile  = open('../ref/ORF')
  for line in ORFfile:
    line = line.rstrip().rsplit('\t')
    ORFdict[line[0]] = [int(line[1]),int(line[2]),line[3]]
    ORFlist.append(line[0])
  ORFfile.close()
  
  ampdict  = {}
  ampfile  = open('../ref/offset')
  for line in ampfile:
    line = line.rstrip().rsplit('\t')
    ampdict[line[0]] = [int(line[1]),int(line[2])]
    for lib in depdict: 
      depdict[lib][line[0]] = 0
      wtdict[lib][line[0]]  = 0
  ampfile.close()
  
  #clean up ORFdict
  seqrange = []
  for amp in ampdict:
    seqrange.extend(range(ampdict[amp][0],ampdict[amp][1]))
  seqrange = list(set(seqrange))
  for ORF in ORFdict:
    ORFrange = range(ORFdict[ORF][0],ORFdict[ORF][1])
    if set(seqrange).isdisjoint(ORFrange):
      ORFlist.remove(ORF)
  
  #reading data
  sumdict['indels'] = 0
  sumdict['extra barcode'] = 0
  sumdict['multiple mutations'] = 0
  sumdict['effective reads'] = 0
  infiles = sorted(glob.glob(workpath+'*.m'))
  #infiles = [workpath+'316.m']
  for infile in infiles:
    print infile
    handle = open(infile)
    for line in handle:
      line = line.rstrip().rsplit('\t')
      if len(line) < 3:
        line = ''.join(line)
        line = line.rsplit(': ')
        if line[0] not in sumdict: sumdict[line[0]] = 0
        sumdict[line[0]] += int(line[1])
        continue
      amp = line[0] ; bc = line[1]; muts = line[2]
      if bc not in libdict: 
        sumdict['extra barcode'] += 1
        continue
      lib = libdict[bc]
      #if lib == 'P62' and amp == 'MF10-A1': print muts
      depdict[lib][amp] += 1
      if muts == 'WT': 
        wtdict[lib][amp] += 1
        sumdict['effective reads'] += 1
        continue
      if '-' in muts:
        sumdict['multiple mutations'] += 1
      muts = muts.rsplit('-')
      if len(muts) > 4: 
        sumdict['indels'] += 1
        depdict[lib][amp] -= 1
        continue
      for mut in muts:
        mut = anchor_mut(mut,amp,ampdict)
        if mut not in genodict: genodict[mut] = {}
        if lib not in genodict[mut]: genodict[mut][lib] = 0
        genodict[mut][lib] += 1
      sumdict['effective reads'] += 1
    handle.close()
  
  #print depdict['P62']['MF10-A1']
  #writefile
  depfile.write('Library\tAmplicon\tWild type\tTotal depth\n')
  for lib in depdict:
    for amp in depdict[lib]:
      depfile.write(lib+'\t'+amp+'\t'+str(wtdict[lib][amp])+'\t'+str(depdict[lib][amp])+'\n')
  depfile.close()
  
  for name in sumdict:
    sumfile.write(name+': '+str(sumdict[name])+'\n')
  sumfile.close()
  
  outfile.write('Pos\tGenotype')
  for ORF in ORFlist:
    outfile.write('\t'+ORF)
  for lib in liblist:
    outfile.write('\t'+lib+'\t'+lib+'_wt\t'+lib+'_depth')
  outfile.write('\n')
  for genotype in genodict:
    pos = genotype[1:-1]
    outfile.write(pos+'\t'+genotype)
    pos = int(pos)
    for ORF in ORFlist:
      ORFrange = range(ORFdict[ORF][0],ORFdict[ORF][1])
      if pos not in ORFrange: 
        outfile.write('\tX')
      else:
        start_pos = ORFrange[0]
        aa_geno   = genotype[0]+str(pos-start_pos)+genotype[-1]
        aamut = translate_mut(aa_geno,ORFdict[ORF][2])
        outfile.write('\t'+aamut)
    for lib in liblist:
      count_wt    = 0
      count_total = 0
      for amp in ampdict:
        if pos <= ampdict[amp][1]+1 and pos >= ampdict[amp][0]+1:
          count_wt    += wtdict[lib][amp]
          count_total += depdict[lib][amp]
      #if '8031' in genotype and lib == 'P62' : print count_wt 
      if lib not in genodict[genotype]: genodict[genotype][lib] = 0
      outfile.write('\t'+str(genodict[genotype][lib])+'\t'+str(count_wt)+'\t'+str(count_total))
    outfile.write('\n')
  outfile.close()
      
if __name__ == '__main__':
  main()
