#!/usr/bin/env python3

import argparse
import gzip
from collections import OrderedDict
import numpy as np
import os
import pandas as pd
# import pyfaidx
import re
import sys

## -----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-k', action='store', dest='kmers', default=None,
                    help='Input file containing k-mers to count/locate.' +
                         ' Either tsv file with kmer column or simple list of k-mers,' +
                         ' one per line. Alternately, can provide comma-delimited'
                         ' (no spaces) of k-mers directly on command line.')
parser.add_argument('-c', action='store', dest='clusters', default=None,
                    help='Alternate input file containing k-mer clusters' +
                         ' to count/locate. Only one of -k or -c should be provided.')
parser.add_argument('-f', action='store', dest='fasta', default=None,
                    help='Input fasta file containing sequences in which k-mers' +
                         ' are to be counted or located.')
parser.add_argument('-d', action='store_true', dest='directional', default=False,
                    help='(d)irectional mode: count k-mers only in the' +
                         ' forward (i.e., not reverse-complemented) orientation.')
parser.add_argument('-l', action='store_true', dest='locate', default=False,
                    help='(l)ocate mode: return tabular output with columns for' +
                         'seqid, location, and kmer/cluster.')
parser.add_argument('-o', action='store_true', dest='overlap', default=False,
                    help='(o)verlapping mode: count overlapping matches as' +
                         'separate instances (not available in locate mode).')

parsed = parser.parse_args()

## -----------------------------------------------------------------
def revComp(s):
    if isinstance(s, pd.Series):
        return pd.Series([revComp(sel) for sel in s], index=s.index)
    comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'N' : 'N', 'T' : 'A',
            'a' : 't', 'c' : 'g', 'g' : 'c', 'n' : 'n', 't' : 'a',
            '[' : ']', ']' : '[', '-' : '-'}
    return ''.join([comp[s[i]] for i in range((len(s)-1), -1, -1)])

def fastaToSeries(infasta):
    # infasta = pyfaidx.Fasta(infasta)
    # return pd.Series({i : str(infasta[i]) for i in infasta.keys()})
    lines = None
    if re.match(r'.*\.[gG][zZ]$', infasta):
        with gzip.open(infasta, 'r') as inf:
            lines = [line.decode('utf-8').strip() for line in inf]
    else:
        with open(infasta, 'r') as inf:
            lines = [line.strip() for line in inf]
    lines = pd.Series(lines).str.replace(r'^>(.*)$', r'>>>\1>>>')
    slurped = ''.join(lines.values)
    slurped = re.sub(r'^>>>', '', slurped)
    lines = pd.Series(slurped.split('>>>'))
    return pd.Series(
        lines.iloc[range(1, len(lines), 2)].values,
        index = lines.iloc[range(0, len(lines), 2)].values
    )

def regexCounts(regex, seqs, overlap=False):
    if overlap:
        regex = r'(?=' + regex + ')'
    counts = seqs.str.upper()\
                 .str.replace(regex, '_')\
                 .str.replace(r'[^_]+', '')\
                 .str.len()
    return counts

def kmerCounts(kmer, seqs, directional=False, overlap=False):
    if not directional:
        return regexCounts(kmer + '|' + revComp(kmer), seqs, overlap)
    else:
        return regexCounts(kmer, seqs, overlap)

def clusterCounts(kmers, seqs, directional=False, overlap=False):
    if not directional:
        kmers = sorted(set(list(kmers) + [revComp(km) for km in kmers]))
    pattern = '|'.join(kmers)
    return regexCounts(pattern, seqs, overlap)

## -----------------------------------------------------------------
def regexLocate(regex, seqs):
    inseqs = seqs
    seqs = '_' + seqs.copy()
    frags = seqs.str.upper()\
                .str.replace('('+regex+')', r'_\1')\
                .str.split('_', expand=True)\
                .stack()\
                .reset_index()
    frags.columns = ['seqid', 'frag_id', 'frag_len']
    frags['frag_len'] = frags['frag_len'].str.len()
    frags.sort_values(['seqid', 'frag_id'], inplace=True)
    frags['location'] = frags[['seqid', 'frag_len']].groupby('seqid')\
                                                    .cumsum()\
                                                    .iloc[:, 0]
    # offset = frags[['seqid', 'location']].groupby('seqid').agg(np.min).iloc[:, 0] -\
    #          frags[['seqid', 'frag_len']].groupby('seqid').agg(np.min).iloc[:, 0]
    # frags['location'] = frags['location'] -\
    #                     offset.reindex(frags['seqid'].values).values
    frags = frags.loc[(frags['frag_id'] > 0) &
                      (frags['location'] <
                       inseqs.reindex(frags['seqid'].values).str.len().values)]
    return pd.DataFrame({
        'seqid' : frags['seqid'],
        'matchid' : frags['frag_id']-1,
        'location' : frags['location']
    })[['seqid', 'matchid', 'location']]

def locateKmers(kmers, seqs, directional=False):
    patterns = {km : km for km in kmers}
    if not directional:
        patterns = {km : km+'|'+revComp(km) for km in kmers}
    locations = {km : regexLocate(patterns[km], seqs)[['seqid', 'location']]
                 for km in kmers}
    for km in kmers:
        locations[km]['kmer'] = km
    return pd.concat(locations).sort_values(['seqid', 'location', 'kmer'])

def locateClusters(clusters, seqs, directional=False):
    out = []
    for cluster in clusters:
        kmers = clusters[cluster]
        if not directional:
            kmers = sorted(set(list(kmers) + [revComp(km) for km in kmers]))
        pattern = '|'.join(kmers)
        locations = regexLocate(pattern, seqs)[['seqid', 'location']]
        locations['cluster'] = cluster
        out.append(locations)
    return pd.concat(out).sort_values(['seqid', 'location', 'cluster'])

## -----------------------------------------------------------------
if parsed.clusters is not None:
    if os.path.exists(parsed.clusters):
        with open(parsed.clusters, 'r') as clusthandle:
            clustLines = [line.strip() for line in clusthandle]
        clusters = OrderedDict()
        activeCluster = None
        for line in clustLines:
            if line.startswith('Cluster '):
                activeCluster = re.sub(r'^Cluster \d+:\s+', '', line) + '+'
                clusters[activeCluster] = []
            elif line != '':
                clusters[activeCluster].append(re.sub(r'-|\s+', '', line))
    else:
        clusters = {parsed.clusters : np.unique(parsed.clusters.split(','))}
else:
    if os.path.exists(parsed.kmers):
        kmers = pd.read_csv(parsed.kmers, sep='\t', header=0, index_col=None)
        if not 'kmer' in kmers.columns:
            with open(parsed.kmers, 'r') as kmhandle:
                kmers = np.unique([line.strip() for line in kmhandle])
        else:
            kmers = kmers['kmer'].unique()
    else:
        kmers = np.unique(parsed.kmers.split(','))
    if not parsed.directional:
        kmers = np.unique([min(km, revComp(km)) for km in kmers])

seqs = fastaToSeries(parsed.fasta)

try:
    if parsed.locate:
        if parsed.clusters is not None:
            out = locateClusters(clusters, seqs, parsed.directional)
        else:
            out = locateKmers(kmers, seqs, parsed.directional)
        out.to_csv(sys.stdout, sep='\t', header=True, index=False)
    elif parsed.clusters is not None:
        pd.DataFrame({
            cl : clusterCounts(clusters[cl], seqs,
                               parsed.directional, parsed.overlap)
            for cl in clusters
        }).reset_index().to_csv(sys.stdout, sep='\t', header=True, index=False)        
    else:
        pd.DataFrame({
            km : kmerCounts(km, seqs,
                            parsed.directional, parsed.overlap)
            for km in kmers
        }).reset_index().to_csv(sys.stdout, sep='\t', header=True, index=False)
except IOError:
    try:
        sys.stdout.close()
    except IOError:
        pass
    try:
        sys.stderr.close()
    except IOError:
        pass
