#!/usr/bin/env python3

import argparse
import itertools
import numpy as np
import os
import pandas as pd
import re
import tempfile

from sarks import Sarks


## -----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-c', action='store', dest='catfasta', default=None)
parser.add_argument('-f', action='store', dest='fasta', default=None)
parser.add_argument('-g', action='store', dest='ginis', default='1.1')
parser.add_argument('-n', action='store_true', dest='negate', default=False)
parser.add_argument('-o', action='store', dest='out', default=None)
parser.add_argument('-r', action='store', dest='reps',
                    type=int, default=100)
parser.add_argument('-p', action='store_true', dest='prune', default=False)
parser.add_argument('-s', action='store', dest='scores', default=None)
parser.add_argument('-z', action='store', dest='nsigma',
                    type=float, default=4.0)
parser.add_argument('-w', action='store', dest='windows',
                    default='250,500,1000,2500')
parser.add_argument('-l', action='store', dest='spatials',
                    default='0')
parser.add_argument('-e', action='store', dest='seed', default=None)

parsed = parser.parse_args()

catfasta = parsed.catfasta
if catfasta is None:
    catfastaObj = tempfile.NamedTemporaryFile()
    catfasta = catfastaObj.name
   
outdir = re.sub(r'/$', '', parsed.out) + '/'

minGinis = [float(g) for g in parsed.ginis.split(',')]
halfWindows = [int(w) for w in parsed.windows.split(',')]
spatialLengths = [int(x) for x in parsed.spatials.split(',')]

theSeed = parsed.seed
if theSeed is not None:
    theSeed = int(theSeed)


## -----------------------------------------------------------------
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

peaks = sarks.multiWindowPeaks(permDists['theta'],
                               prune=parsed.prune, extend=parsed.prune)
for pk in peaks:
    for duple in pk:
        peaks[pk][duple[0]] = duple[1]

peaks = pd.concat(peaks, ignore_index=True)
peaks.to_csv(outdir + 'peaks.tsv', sep='\t', index=False, header=True)


## -----------------------------------------------------------------
## merge spatial peaks
mergedPeaks = []
for i in range(permDists['theta'].shape[0]):
    filt = permDists['theta'].iloc[i]
    sarks = Sarks(parsed.fasta,
                  catFasta = catfasta,
                  scores = scores,
                  halfWindow = int(filt['halfWindow']),
                  regenerateCatFasta = True)
    pks = sarks.peaks(theta = filt['theta'],
                      spatialLength = int(filt['spatialLength']),
                      spatialTheta = filt['spatialTheta'],
                      minGini = filt['minGini'],
                      minSpatialGini = filt['minSpatialGini'])
    if len(pks) > 0:
        pks = sarks.subtable(pks)
        pks = pks.loc[~pks['kmer'].str.match(r'^.*\$')]
        mkm = None
        if filt['spatialLength'] > 0:
            subpeaks = sarks.spatialSubPeaks(pks,
                                             theta = filt['spatialTheta'],
                                             minGini = filt['minGini'])
            if subpeaks.shape[0] > 0:
                mergedKmers = Sarks.mergeKmers(subpeaks)
                if mergedKmers.shape[0] > 0:
                    mkm = pd.DataFrame({'kmer' : sorted(set(mergedKmers['kmer']))})
        else:
            mkm = pd.DataFrame({'kmer' : sorted(set(pks['kmer']))})
        if mkm is not None:
            for idx in filt.index:
                mkm[idx] = filt.loc[idx]
            mergedPeaks.append(mkm.loc[:, list(filt.index) + ['kmer']])

if len(mergedPeaks) > 0 and np.any(np.array(spatialLengths) > 0):
    mergedPeaks = pd.concat(mergedPeaks)
    mergedPeaks.to_csv(outdir + 'merged_peaks.tsv',
                       sep='\t', index=False, header=True)


## -----------------------------------------------------------------
if parsed.catfasta is None:
    catfastaObj.close()
