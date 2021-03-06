#!/usr/bin/env python
def main():
  infile = open('result/readcount.txt')
  outfile = open('result/frequency.txt','w')
  header = infile.readline().rstrip().rsplit('\t')
  outfile.write('\t'.join(header[:3])+'\t'+'\t'.join(header[3::3])+'\n')
  for line in infile:
    line = line.rstrip().rsplit('\t')
    outfile.write('\t'.join(line[:3]))
    for i, count in enumerate(line[3::3]):
      freq = float(count)/float(line[3*i+5])
      outfile.write('\t'+str(freq))
    outfile.write('\n')
  infile.close()
  outfile.close()

if __name__ == '__main__':
  main()
