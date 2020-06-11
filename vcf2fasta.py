#!/usr/bin/env python

import os
import sys
import gzip
import argparse

def get_genotype(genotype, ref, alt):
  out = []
  for g in genotype:
    tmp = g.split(':')[0]
    if tmp == "0/0":
      out.append(''.join((ref,ref)))
    if tmp == "0|0":
      out.append(''.join((ref,ref)))
    elif tmp == "0/1":
      out.append(''.join((ref,alt)))
    elif tmp == "0|1":
      out.append(''.join((ref,alt)))
    elif tmp == "1/1":
      out.append(''.join((alt,alt)))
    elif tmp == "1|1":
      out.append(''.join((alt,alt)))
    elif tmp == ("./."):
      out.append(''.join(('N','N')))
    elif tmp == (".|."):
      out.append(''.join(('N','N')))
  return ';'.join(out)


def main(vcffile):
  N=0
  NBI=0
  NMU=0
  GEN=[]

  with gzip.open(vcffile, 'rb') as f:
    for line in f:
      line=bytes.decode(line)
      if line.startswith('#CHROM'):
        strains = line.split('\t')[9:]
      elif not line.startswith('##'):
        N = N + 1
        REF = line.split('\t')[3]
        ALT = line.split('\t')[4]
        if len(ALT)==1 and ALT != "*":
          NBI = NBI + 1
          GEN.append(get_genotype(line.split('\t')[9:], REF, ALT))
        if len(ALT)>1:
          NMU = NMU + 1
        if N % 1000000 == 0:
          sys.stderr.write('==> Processed %i,000,000 sites\n' % (N/1000000))
        elif N % 100000 == 0:
          sys.stderr.write('==> Processed %i00,000 sites\n' % (N/100000))
        elif N % 10000 == 0:
          sys.stderr.write('==> Processed %i0,000 sites\n' % (N/10000))

  sys.stderr.write('\n==> Found %i bi-allelic SNPs.\n' % (NBI))
  sys.stderr.write('==> Found %i multi-allelic SNPs. These are ignored.\n' % (NMU))

  sys.stderr.write('==> Writing sequences to FASTA file %s\n\n' % os.path.abspath(str(vcffile) + '.fa'))
  with open(str(vcffile) + '.fa', 'a') as f2:
    for i in range(0,len(strains)):
      f2.write('>'+strains[i]+'\n'+''.join([item.split(';')[i] for item in GEN])+'\n')
  f2.close()
      
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='Convert VCF to FASTA file.', \
	  usage = 'vcf2fasta <vcf>')
  parser.add_argument('vcf', help='VCF file', metavar='vcf')
  options = parser.parse_args()
  
  main(options.vcf)
