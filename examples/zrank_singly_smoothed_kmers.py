#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import sys

sarksDir = ''
if len(sys.argv) > 1:
    sarksDir = sys.argv[1]

if len(sarksDir) > 0 and os.path.exists(sarksDir):
    peaksFile = sarksDir + '/peaks.tsv'
    if os.path.exists(peaksFile+'.gz') and not os.path.exists(peaksFile):
        peaksFile += '.gz'
    permFile = sarksDir + '/permdists_windowed.tsv'
    if os.path.exists(permFile+'.gz') and not os.path.exists(permFile):
        permFile += '.gz'

peaks = pd.read_table(peaksFile, header=0, index_col=None)
peaks = peaks.loc[peaks['spatialLength'] == 0]

perms = pd.read_table(permFile, header=0, index_col=None)
perms = perms.loc[perms['spatialLength'] == 0]
permMeans = perms[['halfWindow', 'minGini', '1.0']].groupby(['halfWindow', 'minGini'])\
                                                   .agg(np.mean).iloc[:, 0]
permSds = perms[['halfWindow', 'minGini', '1.0']].groupby(['halfWindow', 'minGini'])\
                                                 .agg(np.std).iloc[:, 0]

peaks['zmax'] = (peaks['windowed'] -
                 permMeans.loc[list(zip(peaks['halfWindow'], peaks['minGini']))].values) /\
                permSds.loc[list(zip(peaks['halfWindow'], peaks['minGini']))].values

peaks.sort_values('zmax', ascending=False, inplace=True)

deduped = peaks[['kmer', 'halfWindow', 'minGini', 'zmax']].drop_duplicates('kmer')
deduped['halfWindow'] = deduped['halfWindow'].astype(int)
deduped.to_csv(sys.stdout, sep='\t', header=True, index=False)