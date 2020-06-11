#!/usr/bin/env python


#%% MODULES TO IMPORT 
from collections import defaultdict
import os
import sys
import gzip
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#%% GET SAMPLE NAMES FROM VCF FILE        
def get_names(vcffile):
  with gzip.open(vcffile, 'rb') as f:
    for line in f:
      line=bytes.decode(line)
      if line.startswith('#CHROM'):
        line = line.strip()
        return line.split('\t')
        break


#%% PLOT HISTOGRAM OF GENOME-WIDE FREQUENCIES        
def plot_histogram(frequencies, name):
  n, bins, patches = plt.hist(frequencies, 50, facecolor='g', alpha=0.75)
  plt.xlabel('Alternate allele frequency')
  plt.ylabel('Probability')
  plt.title(str(name))
  plt.grid(False)
  plt.savefig('FreqHistograms/' + str(name) + '.hist.pdf')
  plt.clf()
  plt.close()


#%% PLOT HISTOGRAM AND SCATTERPLOTS OF FREQUENCIES PER CHROMOMSOME
def plot_scatterplot(freqA, freqB, name, chromosome):
  xall = list(map(int, freqA.keys())) + list(map(int, freqB.keys()))
  yall = list(freqA.values()) + list(freqB.values())
  y = freqA.values() 
  
  plt.figure(figsize=(20,6))
  plt.title(str(name)+ '_' + chromosome)
  
  gs = gridspec.GridSpec(1, 2, width_ratios=[5,1])
  ax_main = plt.subplot(gs[0])
  ax_main.set_ylim([-0.05,1.05])
  ax_yDist = plt.subplot(gs[1], sharey=ax_main)
  ax_yDist.set_ylim([-0.05,1.05])

  try:
    ax_main.set_xlim([0, np.max(xall)])
    ax_main.scatter(xall, yall)
    ax_main.set(xlabel="index", ylabel="Alternate allele frequency")
    ax_yDist.hist(y, bins=50, orientation='horizontal', facecolor='g', align='mid')
    ax_yDist.set(xlabel='Frequency')
  except ValueError:
    sys.stderr.write('==> No heterozygous sites in chromosome %s of isolate %s\n' % (chromosome, name))

  plt.savefig('FreqScatterplots/' + str(name) + '/' + str(name) + '_' + chromosome +'.scatter.pdf')
  plt.clf()
  plt.close()



#%% RUNS MAIN MODULE
#vcffile='/Users/fvandenbroeck/Documents/6_Science/1_Projects/2018/4_Lperuviana/3_SNPsINDELs/1_VCFfiles/LBRACOMPLEX.GENO.SNP.1500.ANN.PASS.g.vcf.gz'
#sample='HB22A1'
def main(vcffile, sample):
  N=0
  NBI=0
  NMU=0
  SNPCOLUMN = get_names(vcffile).index(sample)
  chromosomes=[]
  HETfreqs=defaultdict(lambda: defaultdict(list))
  HOMfreqs=defaultdict(lambda: defaultdict(list))

  sys.stderr.write('==> Processing sites\n')
  with gzip.open(vcffile, 'rb') as f:
    for line in f:
      line=bytes.decode(line)
      if not line.startswith(('#','LbrM_00','LbrM_maxicircle')):
        N = N + 1
        CHR = line.split('\t')[0]
        POS = line.split('\t')[1]
        ALT = line.split('\t')[4]
        FILTER = line.split('\t')[6]
        SNP = line.split('\t')[SNPCOLUMN]

        if CHR not in chromosomes:
          chromosomes.append(CHR)
        
        if len(ALT)==1 and FILTER == 'PASS':     # keeps on bi-allelic SNPs
          NBI = NBI + 1
          frq = SNP.split(':')
          if frq[0] == '0/1':       # only alt freq of heterozygous sites
            try:
                HETfreqs[CHR][POS] = (float(frq[1].split(',')[1])/float(frq[2]))
            except ZeroDivisionError:
                sys.stderr.write('WARNING: %s has zero depth at SNP located on chromosome %s and position %s\n' % (sample, CHR, POS))

          if frq[0] in ('0/0', '1/1'):       # only alt freq of homozygous sites
            try:
                HOMfreqs[CHR][POS] = (float(frq[1].split(',')[1])/float(frq[2]))
            except ZeroDivisionError:
                sys.stderr.write('WARNING: %s has zero depth at SNP located on chromosome %s and position %s\n' % (sample, CHR, POS))

        if len(ALT)>1:
          NMU = NMU + 1

  sys.stderr.write('==> Found %i bi-allelic SNPs.\n' % (NBI))
  sys.stderr.write('==> Found %i multi-allelic SNPs. These are ignored.\n' % (NMU))


#%% CREATES HISTOGRAMS AND SCATTERPLOTS
  sys.stderr.write('==> Writing plots to disk...\n')
  
  if not os.path.exists('FreqHistograms'):
    os.makedirs('FreqHistograms')  
  if not os.path.exists('FreqScatterplots'):
    os.makedirs('FreqScatterplots')
  if not os.path.exists('FreqScatterplots/' + str(sample) + '/'):
    os.makedirs('FreqScatterplots/' + str(sample) + '/')
  
  gwfreqs = []
  for c in range(0, len(chromosomes), 1):
    plot_scatterplot(freqA = HETfreqs[chromosomes[c]], freqB = HOMfreqs[chromosomes[c]], name = sample, chromosome = chromosomes[c])
    gwfreqs += HETfreqs[chromosomes[c]].values()
#    for r in HETfreqs[chromosomes[c]].values():
#      gwfreqs.append(r)
  plot_histogram(frequencies = gwfreqs, name = sample)


#%% MAIN
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='Creates histograms and scatterplots of alternate allele read depth frequencies at heterozygous sites.', \
      usage = 'vcf2freq.py <name> <vcf>')
  parser.add_argument('name', help='Sample name for which frequencies should be plotted', metavar='vcf')
  parser.add_argument('vcf', help='VCF file', metavar='vcf')
  options = parser.parse_args()
  
  main(options.vcf, options.name)