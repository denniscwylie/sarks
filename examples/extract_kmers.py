#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

## -----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-d', action='store_true', dest='directional', default=False,
                    help='(d)irectional mode: do not consider k-mer equivalent' +
                         ' to its reverse-complement')
parser.add_argument(action='store', dest='input', default=None,
                    help='tsv file with kmer column')

parsed = parser.parse_args()

## -----------------------------------------------------------------
def revComp(s):
    if isinstance(s, pd.Series):
        return pd.Series([revComp(sel) for sel in s], index=s.index)
    comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'N' : 'N', 'T' : 'A',
            '[' : ']', ']' : '[', '-' : '-'}
    return ''.join([comp[s[i]] for i in range((len(s)-1), -1, -1)])

## -----------------------------------------------------------------
if parsed.input is not None:
    tab = pd.read_csv(parsed.input, sep='\t', header=0, index_col=None)
    kmers = tab['kmer'].unique()
    if not parsed.directional:
        kmers = np.unique([min(km, revComp(km)) for km in kmers])
    for km in kmers:
        print(km)
