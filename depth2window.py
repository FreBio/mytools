#!/usr/bin/env python

'''
Input = samtools depth file
Output = median depth per 2kb windows
'''

from __future__ import division
from collections import defaultdict

import os
import sys
import gzip
import argparse
import numpy as np


def main(outfile, infile, length):

  sys.stderr.write('\n==> Reading and storing depth file. This might take a while.')
  depths = defaultdict(list)
  positions = defaultdict(list)
  with gzip.open(infile, 'rb') as f:
    for line in f:
      line=bytes.decode(line)
      if line.startswith(('Tb927_0','Tb927_10_v5.1','Tb927_11_v5.1')):
        loc=line.split()
        positions[str(loc[0])].append(int(loc[1]))
        depths[str(loc[0])].append(int(loc[2])) 
 
  sys.stderr.write('\n==> Estimating median read depths in %s bp windows.' % (length))
  
  chromosomes = depths.keys()  
  sorted(chromosomes)
  
  with open(outfile + '_chr_depths.txt', 'a') as fd:
    fd.write('%s %s\n' % ('CHR', 'MEDIANCHRDEPTH'))
    for chr in chromosomes:
      fd.write('%s %s\n' % (chr, np.median(np.array(depths[chr]))))
    fd.close()

  with open(outfile + '_window_depths.txt', 'a') as f:
    f.write('%s %s %s %s %s %s %s\n' % ('CHR', 'BIN', 'N_ZERO_DEPTH', 'MEDIANHAPDEPTH', 'MEANHAPDEPTH', 'MEDIANDEPTH', 'MEANDEPTH'))
    for chr in chromosomes:
      chrbins = np.arange(0, int(positions[chr][-1]), int(length))
      chrpos = np.array(positions[chr])
      chrdepth = np.array(depths[chr])
      chrdepthmed = np.median(chrdepth)

      for i in range(0, len(chrbins)-1):
        depthsbin = chrdepth[np.logical_and(chrpos > chrbins[i], chrpos <= chrbins[i+1])]
        f.write('%s %s %s %s %s %s %s\n' % (chr, i*int(length), len(depthsbin[depthsbin==0]), round(np.median(depthsbin)/chrdepthmed, 2), round(np.mean(depthsbin)/chrdepthmed, 2), np.median(depthsbin), np.mean(depthsbin)))
    f.close()
    
  sys.stderr.write('\n==> Results written to %s \n\n' % (os.path.abspath(outfile + '_windowdepth.txt')))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='Median depth in bp windows', \
	  usage = 'depth2kb [options] <out> <depthfile>')
  parser.add_argument('out', help='Prefix used for labeling output files', metavar='out')
  parser.add_argument('depthfile', help='Read depths as outputted from samtools', metavar='depth')
  parser.add_argument('--window', help='Window length in bp [%(default)s].', default=2000, metavar='')
  options = parser.parse_args()
  
  main(options.out, options.depthfile, options.window)
