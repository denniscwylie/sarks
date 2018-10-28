#!/usr/bin/env python3

import argparse
import itertools
import numpy as np
import os
import pandas as pd
import re
from scipy import stats
import tempfile

from sarks import Sarks


## -----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='indir', default=None,
                    help='input directory for sarkstest.py is'+
                         ' output directory from sarkselect.py')
parser.add_argument('-o', action='store', dest='out', default=None,
                    help='output directory to be created/overwritten')
parser.add_argument('-f', action='store', dest='fasta', default=None,
                    help='input fasta file containing sequences to analyze')
parser.add_argument('-s', action='store', dest='scores', default=None,
                    help='input scores tsv file: col1=seqids, col2=numeric scores')
parser.add_argument('-w', action='store', dest='windows',
                    default='250,500,1000,2500',
                    help='half window width (kappa) for first kernel smoothing pass,' +
                         ' can supply multiple values using commas (no spaces)')
parser.add_argument('-l', action='store', dest='spatials',
                    default='0',
                    help='spatial smoothing length (lambda),' +
                         ' can supply multiple values using commas (no spaces)' +
                         ' [default 0]')
parser.add_argument('-g', action='store', dest='ginis', default='1.1',
                    help='parameter for calculation of Gini impurity filter (gamma),' +
                         ' can supply multiple values using commas (no spaces)' +
                         ' [default 1.1]')
parser.add_argument('-r', action='store', dest='reps',
                    type=int, default=100,
                    help='number R_2 of permutations used to test significance thresholds' +
                         ' [default 100]')
parser.add_argument('-z', action='store', dest='nsigma',
                    type=float, default=4.0,
                    help='multiple z of standard deviations above mean' +
                         ' (of maximum smoothed suffix scores obtained after' +
                         ' randomly permuting scores assigned to sequences)' +
                         ' defining threshold (Section S2.6, Eq (S24-S25) of paper)' +
                         ' [default 4.0]')
parser.add_argument('-p', action='store_true', dest='prune', default=False,
                    help='(optional flag) remove redundant kmer output'+
                         ' (Section S2.3 of paper), most useful for small data sets')
parser.add_argument('-n', action='store_true', dest='negate', default=False,
                    help='(optional flag) find motifs associated with' +
                         ' lowest/most negative scores' +
                         ' instead of highest/most positive scores')
parser.add_argument('-c', action='store', dest='catfasta', default=None,
                    help='(optional) path for concatenated fasta file')
parser.add_argument('-e', action='store', dest='seed', default=None,
                    help='(optional) seed for random number generator')

parsed = parser.parse_args()

catfasta = parsed.catfasta
if catfasta is None:
    catfastaObj = tempfile.NamedTemporaryFile()
    catfasta = catfastaObj.name

outdir = re.sub(r'/$', '', parsed.outdir) + '/'
if parsed.indir is not None:
    indir = re.sub(r'/$', '', parsed.indir) + '/'

minGinis = [float(g) for g in parsed.ginis.split(',')]
halfWindows = [int(w) for w in parsed.windows.split(',')]
spatialLengths = [int(x) for x in parsed.spatials.split(',')]

theSeed = parsed.seed
if theSeed is not None:
    theSeed = int(theSeed)


## -----------------------------------------------------------------
if parsed.scores is not None and parsed.fasta is not None:
    scores = pd.read_csv(parsed.scores, sep='\t',
                         index_col=0, header=0).iloc[:, 0]
    if parsed.negate is not None and parsed.negate:
        scores = -scores
    
    sarks = Sarks(parsed.fasta, catfasta, scores,
                  halfWindow = halfWindows[0],
                  regenerateCatFasta = parsed.catfasta is None)
    
    filters = pd.DataFrame(list(itertools.product(
        spatialLengths,
        minGinis
    )), columns=['spatialLength', 'minGini'])
    filters['minSpatialGini'] = filters['minGini']
    
    permDists = sarks.multiWindowPermute(
        halfWindows = halfWindows,
        filters = filters,
        reps = parsed.reps,
        seed = theSeed,
        nsigma = parsed.nsigma
    )
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    for key in permDists:
        permDists[key].to_csv(outdir + 'permdists_' + key + '.tsv',
                              sep='\t', index=False, header=True)
elif os.path.exists(parsed.outdir):
    permDists = {}
    permDists['windowed'] = pd.read_table(
        outdir + 'permdists_windowed.tsv', index_col=None, header=0
    )
    permDists['windowed'][1.0] = permDists['windowed']['1.0']
    permDists['spatial_windowed'] = pd.read_table(
        outdir + 'permdists_spatial_windowed.tsv', index_col=None, header=0
    )
    permDists['spatial_windowed'][1.0] = permDists['spatial_windowed']['1.0']

    
## -----------------------------------------------------------------
if indir is not None:
    inPerm = pd.read_table(indir + 'permdists_windowed.tsv',
                           index_col=None, header=0)
    inSpat = pd.read_table(indir + 'permdists_spatial_windowed.tsv',
                           index_col=None, header=0)
    
    def calcThresh(x):
        return np.mean(x) + parsed.nsigma * np.std(x, ddof=1)
    
    aggCols = ['halfWindow', 'minGini']
    thetas = inPerm[aggCols + ['1.0']].loc[inPerm['spatialLength'] == 0]\
                                      .groupby(aggCols)\
                                      .agg(calcThresh)\
                                      .iloc[:, 0]
    if inSpat.shape[0] > 0:
        aggCols = ['halfWindow', 'spatialLength', 'minGini']
        spatThetas = inSpat[aggCols + ['1.0']].loc[inSpat['spatialLength'] > 0]\
                                              .groupby(aggCols)\
                                              .agg(calcThresh)\
                                              .iloc[:, 0]
    
    singlySmoothed = permDists['windowed'].copy()
    singlySmoothed = singlySmoothed.loc[singlySmoothed['spatialLength'] == 0]
    singlySmoothed['theta'] = thetas.loc[
        list(zip(singlySmoothed['halfWindow'], singlySmoothed['minGini']))
    ].values
    
    posReps = set(singlySmoothed.loc[singlySmoothed[1.0] >=
                                     singlySmoothed['theta'],
                                     'rep'])

    if inSpat.shape[0] > 0:
        spatSmoothed = permDists['spatial_windowed'].copy()
        spatSmoothed = spatSmoothed.loc[spatSmoothed['spatialLength'] > 0]
        spatSmoothed['spatialTheta'] = spatThetas.loc[
            list(zip(spatSmoothed['halfWindow'],
                     spatSmoothed['spatialLength'],
                     spatSmoothed['minGini']))
        ].values
        
        posReps |= set(spatSmoothed.loc[spatSmoothed[1.0] >=
                                        spatSmoothed['spatialTheta'],
                                        'rep'])
    
    def binomCI_exact(x, n, confLevel=0.95):
        lower, upper = 0, 1
        halfAlpha = (1 - confLevel) / 2.0
        if x > 0:
            lower = stats.beta.isf(1-halfAlpha, x, n+1-x)
        if x < n:
            upper = stats.beta.isf(halfAlpha, x+1, n-x)
        return lower, upper

    nTestReps = len(set(permDists['windowed']['rep']) |
                    set(permDists['spatial_windowed']['rep']))
    ci = binomCI_exact(len(posReps), nTestReps)
    print(str(len(posReps)) + ' / ' + str(nTestReps) +
          '\npoint estimate: ' + str(round(100.*len(posReps)/nTestReps, 3)) + '%' +
          '\n95% CI: (' +
          str(round(100.*ci[0], 3)) + '%, ' + str(round(100.*ci[1], 3)) + '%)')
