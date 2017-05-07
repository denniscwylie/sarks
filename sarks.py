#!/usr/bin/env python3

from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import intervaltree
import numpy as np
import os
import pandas as pd
from pyfaidx import Fasta
import re
from scipy import stats
from subprocess import check_output
import tempfile

from windginiimp import windGiniImpurities


class Sarks(object):
    
    def __init__(self, inFasta, catFasta, scores,
                 suffixArrayFile=None,
                 halfWindow=250,
                 minScore=-10, maxScore = 10, maxRun=np.inf,
                 seqs=None, catSeq=None, bounds=None, sa=None,
                 windGini=None, spatialLength=None, spatGini=None,
                 regenerateSuffixArray=False):
        self.inFasta = inFasta
        self.catFasta = catFasta
        self.suffixArrayFile = suffixArrayFile
        self.halfWindow = halfWindow
        self.minScore = minScore
        self.maxScore = maxScore
        self.maxRun = maxRun
        if seqs is not None:
            self.seqs = seqs
        else:
            self.seqs = Fasta(inFasta)
        self.transcripts = sorted(list(set(scores.index) &
                                       set(self.seqs.keys())))
        self.scores = scores.loc[self.transcripts].copy()
        if catSeq is not None and bounds is not None:
            self.catSeq = catSeq
            self.bounds = bounds
        else:
            self.concatenateSeqs(catFasta, regenerateSuffixArray)
        if sa is not None:
            self.sa = sa
        else:
            self.calcSuffixArray(catFasta, suffixArrayFile,
                                 regenerateSuffixArray)
        if windGini is not None:
            self.windGini = windGini
        else:
            self.windowGini(recalculate=True)
        self.spatialLength = spatialLength
        self.spatGini = spatGini
        self.window(recalculate=True)

    def concatenateSeqs(self, catFasta,
                        regenerateSuffixArray=False):
        """concatenate seqs together and write to catFasta"""
        if regenerateSuffixArray or not os.path.exists(catFasta):
            catSeq = SeqRecord(Seq("".join([str(self.seqs[t]) + "$"
                                            for t in self.transcripts]),
                                   DNAAlphabet()),
                               id = "concatenated_promoters",
                               description = "")
            if os.path.exists(catFasta + ".fai"):
                os.remove(catFasta + ".fai")
            with open(catFasta, "w") as outhandle:
                SeqIO.write([catSeq], outhandle, "fasta")
        ## now read catFasta back in and calculate block boundaries
        self.catSeq =  Fasta(catFasta)['concatenated_promoters']
        self.bounds = np.cumsum([len(self.seqs[t]) + 1
                                 for t in self.transcripts])

    def calcSuffixArray(self, catFasta, suffixArrayFile,
                        regenerateSuffixArray=False):
        """run Burrows-Wheeler algorithm to get suffix array"""
        if suffixArrayFile is None:
            suffixArrayFile = re.sub(r"\..*", "", catFasta) + "_sa.txt"
        if regenerateSuffixArray or not os.path.exists(suffixArrayFile + ".gz"):
            check_output(
                "suffix-array " + catFasta +
                " | grep -v '>' > " + suffixArrayFile,
                shell = True
            )
            check_output("gzip -f " + suffixArrayFile, shell=True)
        self.sa = pd.Series(np.loadtxt(suffixArrayFile + ".gz").astype(int))
        ## need to remove extra end-of-string position
        self.sa = self.sa.loc[self.sa < len(self.catSeq)]
        
    def sourceBlock(self, s):
        return np.searchsorted(self.bounds, s+1)

    def sourceScore(self, s, process=False, maxRun=None):
        srcBlk = self.sourceBlock(s)
        out = pd.DataFrame({
            'block' : srcBlk,
            'score' : self.scores.iloc[srcBlk].values
        }, index=s)
        if process:
            if maxRun is None:
                maxRun = self.maxRun
            if maxRun < np.inf:
                out['pos'] = out.index
                out['inrun'] = (out['block'].diff() == 0).astype(int)
                out['rundex'] = np.array(range(out.shape[0])) - \
                                out['inrun'].cumsum().values
                out = out.groupby('rundex').nth(list(range(maxRun)))
                out.index = out['pos']
            out.loc[out['score'] < self.minScore, 'score'] = self.minScore
            out.loc[out['score'] > self.maxScore, 'score'] = self.maxScore
        out['block'] = out['block'].astype(int)
        return out[['block', 'score']]

    def windowGini(self, recalculate=False):
        if "windGini" not in dir(self) or recalculate:
            block = pd.Series(self.sourceBlock(self.sa), index=self.sa)
            self.windGini = windGiniImpurities(block, int(self.halfWindow))
        return self.windGini

    def spatialGini(self, recalculate=False):
        if self.spatGini is None or recalculate:
            sortedGini = self.windGini.copy()
            sortedGini.sort_index(inplace=True)
            self.spatGini = sortedGini.copy()
            giniCumul = self.spatGini.cumsum()
            self.spatGini = giniCumul.iloc[(self.spatialLength-1):] - \
                            np.concatenate([
                                np.zeros(1),
                                giniCumul.iloc[0:(-self.spatialLength)].values
                            ])
            self.spatGini /= self.spatialLength
            self.spatGini.index = \
                    sortedGini.iloc[range(len(self.spatGini))].index
        return self.spatGini
    
    def windowBlockCharEntropy(self, s, j):
        if "__contains__" in dir(j):
            return np.sum([self.windowBlockCharEntropy(s, jel) for jel in j])
        else:
            return stats.entropy(self.kmers(s+j, 1).value_counts(), base=2)

    def window(self, halfWindow=None, maxRun=None, recalculate=False):
        if halfWindow is None:
            halfWindow = self.halfWindow
        if halfWindow != self.halfWindow or recalculate:
            if halfWindow != self.halfWindow:
                self.windowGini(recalculate=True)
            self.halfWindow = halfWindow
            if maxRun is not None:
                self.maxRun = maxRun
            saScores = self.sourceScore(self.sa,
                                        process=True,
                                        maxRun=self.maxRun)['score']
            saCumul = saScores.cumsum()
            windowSize = (2 * halfWindow) + 1
            self.windowed = saCumul.iloc[(windowSize-1):] - \
                            np.concatenate([
                                np.zeros(1),
                                saCumul.iloc[0:(-windowSize)].values
                            ])
            self.windowed /= windowSize
            saWindowCenters = range(
                int(halfWindow),
                int(halfWindow) + len(self.windowed)
            )
            self.windowed.index = saScores.iloc[saWindowCenters].index
        return self.windowed

    def spatialWindow(self, length):
        if 'spatialLength' not in dir(self) or \
           length != self.spatialLength or \
           'spatialWindowed' not in dir(self):
            self.spatGini = None
            self.spatialLength = length
            windowed = self.window().copy()
            windowed.sort_index(inplace=True)
            windCumul = windowed.cumsum()
            self.spatialWindowed = windCumul.iloc[(length-1):] - \
                                   np.concatenate([
                                       np.zeros(1),
                                       windCumul.iloc[0:(-length)].values
                                   ])
            self.spatialWindowed /= length
            self.spatialWindowed.index = \
                    windowed.iloc[range(len(self.spatialWindowed))].index
        return self.spatialWindowed

    def i2s(self, i):
        return int(self.sa.loc[i])

    def s2i(self, s):
        return int(self.sa.loc[self.sa == s].index[0])

    def findKmer(self, kmer):
        if 'kmerBounds' in dir(self):
            if kmer in self.kmerBounds:
                return self.kmerBounds[kmer]
        else:
            self.kmerBounds = {}
        k = len(kmer)
        def km(i):
            return str(self.catSeq[int(self.sa.loc[i]):int(self.sa.loc[i]+k)])
        def startGood(start):
            kms = km(start)
            if start == 0 and kms == kmer:
                return True
            elif kms == kmer and km(start-1) != kmer:
                return True
            else:
                return False
        startLower = 0
        startUpper = len(self.sa)
        start = startLower + (startUpper-startLower) // 2
        while not startGood(start):
            if km(start) >= kmer:
                startUpper = start
            else:
                startLower = start
            newStart = startLower + (startUpper-startLower) // 2
            if newStart == start:
                start = startUpper
            else:
                start = newStart
        def endGood(end):
            kmem1 = km(end-1)
            if end == len(self.sa) and kmem1 == kmer:
                return True
            elif kmem1 == kmer and km(end) != kmer:
                return True
            else:
                return False
        endLower = start
        endUpper = len(self.sa)
        end = endLower + (endUpper-endLower) // 2
        while not endGood(end):
            if km(end) <= kmer:
                endLower = end
            else:
                endUpper = end
            newEnd = endLower + (endUpper-endLower) // 2
            if newEnd == end:
                end = endUpper
            else:
                end = newEnd
        self.kmerBounds[kmer] = (start, end)
        return start, end

    def pruneIntervals(intervals):
        tree = intervaltree.IntervalTree(intervals=[intervaltree.Interval(*iv)
                                                    for iv in intervals])
        nests = tree.find_nested()
        out = set(tree)
        for containing in nests:
            out -= nests[containing]
        return sorted([iv[0:2] for iv in out])

    def peaks(self, theta, prune=True,
              spatialLength=None, spatialTheta=None,
              k=12, minK=None,
              minGini=None, minSpatialGini=None,
              s=None, deduplicate=False):
        scores = None
        if spatialLength is not None:
            scores = self.spatialWindow(spatialLength)
        else:
            scores = self.windowed
        if spatialTheta is None:
            spatialTheta = theta
        scoreDiffB = self.windowed.sort_index().diff()
        scoreDiffF = scoreDiffB.copy().iloc[1:].copy()
        scoreDiffF.index = scoreDiffB.index[:-1]
        scoreDiffB = scoreDiffB.iloc[:-1]
        pos = self.filter(minGini=minGini, minSpatialGini=minSpatialGini,
                          minK=minK,
                          theta=theta,
                          spatialLength=spatialLength, spatialTheta=spatialTheta,
                          nearbyPars=None,
                          k=k, prune=prune, deduplicate=deduplicate, s=s)
        pos = pos.loc[(scoreDiffB.loc[pos] >= 0) &
                      (scoreDiffF.loc[pos] <= 0)]
        return self.sa.loc[self.sa.isin(set(pos))]

    def filter(self, minGini=None, minSpatialGini=None, minK=None,
               theta=None, spatialLength=None, spatialTheta=None,
               nearbyPars=None, k=12, prune=True, deduplicate=True, s=None):
        pos = pd.Series(self.windowed.index)
        pos.index = pos
        posSubtable = None
        if minGini is not None:
            pos = pos.loc[self.windGini[pos] >= minGini]
        if theta is not None:
            pos = pos.loc[(self.windowed.loc[pos] >= theta)]
        if spatialLength is not None:
            spat = self.spatialWindow(spatialLength)
            if minSpatialGini is not None:
                pos = pos.loc[self.spatialGini()[pos] >= minSpatialGini]
            if spatialTheta is not None:
                pos = pos.loc[spat.loc[pos] >= spatialTheta]
        if s is not None:
            pos = pos.loc[pos.isin(s)]
        if deduplicate:
            posSubK = self.quickSubtable(pos, k).copy()
            posSubK.sort_values('windowed', inplace=True)
            posSubK = posSubK.groupby('kmer').first()
            pos = pos.loc[posSubK['s']]
        if prune:
            posSubtable = self.subtable(pos, k).loc[pos]
            posIntervals = [(ess, ess+len(posSubtable.loc[ess, 'kmer']))
                            for ess in posSubtable.index]
            posIntervals = Sarks.pruneIntervals(posIntervals)
            pos = pos.loc[pos.isin([iv[0] for iv in posIntervals])]
            posSubtable = posSubtable.loc[pos]
        if minK is not None:
            if posSubtable is None:
                posSubtable = self.subtable(pos, k).loc[pos]
            pos = pos.loc[posSubtable['khat'] >= minK]
        if nearbyPars is not None:
            if posSubtable is None:
                posSubtable = self.subtable(pos, k).loc[pos]
            pos = pd.Series([s for s in pos if self.rankNearbyKmers(
                posSubtable.loc[s, 'kmer'],
                **nearbyPars
            ).shape[0] > 0])
            pos.index = pos
        return pos

    def extendKmers(self, subtable):
        subtable = subtable.copy()
        maxLen = subtable['kmer'].str.len().max()
        oldS = pd.Series(-1, index=subtable.index)
        while not (oldS == subtable['s']).all():
            oldS = subtable['s']
            kmerSet = set(subtable['kmer'])
            kmerLens = subtable['kmer'].str.len()
            for s in subtable['s'].loc[kmerLens < maxLen].index:
                sKmerLen = len(subtable.loc[s, 'kmer'])
                shiftInt = (s, s+sKmerLen)
                for ext in range(maxLen-sKmerLen, 0, -1):
                    for shift in range(-ext, 1):
                        shiftInt = (int(s)+shift, int(s+sKmerLen+ext+shift))
                        shiftKmer = str(self.catSeq[shiftInt[0]:shiftInt[1]])
                        if shiftKmer in kmerSet:
                            subtable.loc[s, 'i'] += shift
                            subtable.loc[s, 's'] += shift
                            subtable.loc[s, 'wi'] += shift
                            subtable.loc[s, 'kmer'] = shiftKmer
                            subtable.loc[s, 'khat'] = np.nan
                            subtable.loc[s, 'gini'] = np.nan
                            subtable.loc[s, 'windowed'] = np.nan
                            if "spatialWindowed" in dir(self):
                                subtable.loc[s, 'spatial_windowed'] = np.nan
                            break
                    else:
                        continue
                    break
        subtable.index = subtable['s']
        return subtable

    def similarPositions(self, s, d):
        if "__contains__" not in dir(s):
            s = [s]
        if "__contains__" not in dir(d):
            d = (-d, d+1)
        eyes = list(set().union(*[
            set(range(i+d[0], i+d[1]))
            for i in self.sa.loc[self.sa.isin(s)].index
        ]))
        return self.sa.loc[eyes]

    def nearbyPositions(self, s, d):
        if "__contains__" not in dir(s):
            s = [s]
        if "__contains__" not in dir(d):
            d = (0, d)
        esses = list(set().union(*[
            set(range(ess+d[0], ess+d[1]))
            for ess in s
        ]))
        return self.sa.loc[self.sa.isin(esses)]

    def kmers(self, s, k=12, k0=0, sanitize=True):
        catSeqMaxPos = len(self.catSeq) - 1
        def proc1(i):
            return max(0, min(catSeqMaxPos-1, int(i)))
        def proc2(i):
            return max(0, min(catSeqMaxPos, int(i)))
        out = pd.Series({
            w : str(self.catSeq[proc1(w+k0):proc2(w+k)])
            for w in s
        })
        if sanitize:
            out = out.loc[out.str.len() == (k-k0)]
        return out

    def prefixAgreeSum(s, esses):
        esses = pd.Series(esses)
        out = 0
        aliveIndices = esses.index
        for i, letter in enumerate(str(s)):
            iLetters = esses[aliveIndices].str[i].dropna()
            aliveIndices = esses[iLetters.index][iLetters == letter].index
            out += len(aliveIndices)
        return out / len(esses)

    def subtable(self, s, k=12):
        srcScr = self.sourceScore(s)
        eyes = [self.s2i(sel) for sel in s]
        kmers = self.kmers(s, k, sanitize=False).loc[s]
        block = self.scores.iloc[srcScr.loc[s, 'block']].index
        wi = s.values - self.bounds[srcScr.loc[s, 'block'].values]
        wi = wi + np.array([len(self.seqs[t]) + 1 for t in block])
        out = pd.DataFrame({
            "i" : eyes,
            "s" : s.values,
            "kmer" : kmers,
            "khat" : [Sarks.prefixAgreeSum(
                km,
                self.kmers(self.sa[
                    list(range(i-self.halfWindow, i)) +
                    list(range(i+1, i+self.halfWindow+1))
                ], k, 0, False)) for i, km in zip(eyes, kmers)],
            "block" : block,
            "wi" : wi,
            "gini" : self.windGini.loc[s],
            "score" : srcScr.loc[s, 'score'],
            "windowed" : self.windowed.loc[s]
        }, index=s)[['i', 's', 'kmer', 'khat', 'block',
                     'wi', 'gini', 'score', 'windowed']]
        if "spatialWindowed" in dir(self):
            out['spatial_windowed'] = self.spatialWindowed[s]
        for idx in out.index:
            kidx = int(np.round(out.loc[idx, 'khat']))
            out.loc[idx, 'kmer'] = out.loc[idx, 'kmer'][0:kidx]
        return out

    def quickSubtable(self, s, k=12):
        srcScr = self.sourceScore(s)
        eyes = [self.s2i(sel) for sel in s]
        kmers = self.kmers(s, k, sanitize=False).loc[s]
        block = self.scores.iloc[srcScr.loc[s, 'block']].index
        wi = s.values - self.bounds[srcScr.loc[s, 'block'].values]
        wi = wi + np.array([len(self.seqs[t]) + 1 for t in block])
        out = pd.DataFrame({
            'i' : eyes,
            's' : s.values,
            'kmer' : kmers,
            'block' : block,
            'wi' : wi,
            'gini' : self.windGini.loc[s],
            'score' : srcScr.loc[s, 'score'],
            'windowed' : self.windowed.loc[s]
        }, index=s)[['i', 's', 'kmer', 'block',
                     'wi', 'gini', 'score', 'windowed']]
        return out

    def kmerMatrix(self, s, k=12, k0=0):
        kmers = self.kmers(s, k, k0=k0, sanitize=True).loc[s]
        alignMatrix = pd.DataFrame([list(str(km)) for km in kmers])
        return alignMatrix

    def kmerPositionFrequencies(self, s, k=12, k0=0, alphabet=None):
        if alphabet is None:
            alphabet = ['A', 'C', 'G', 'T']
        counts = pd.DataFrame(np.zeros((k-k0, len(alphabet))),
                              columns = alphabet)
        counts = self.kmerMatrix(s, k=k, k0=k0)
        for letter in alphabet:
            counts[letter] = (alignMatrix == letter).sum(axis=0)
        return counts

    def nearbyMatrix(self, seed, theta, k=12, k0=0):
        s = self.sa.loc[range(*self.findKmer(seed))].copy()
        s = s.iloc[(self.scores.iloc[self.sourceBlock(s)] >= theta).values]
        kmers = self.kmers(s, k, k0=k0, sanitize=True).loc[s]
        s = s.iloc[(~kmers.str.match(r".*\$")).values]
        return self.kmerMatrix(s, k, k0=k0)

    def nearbyPositionFrequencies(self, s, k=12, k0=0, alphabet=None):
        if alphabet is None:
            alphabet = ['A', 'C', 'G', 'T']
        counts = self.nearbyMatrix(s, k=k, k0=k0)
        for letter in alphabet:
            counts[letter] = (alignMatrix == letter).sum(axis=0)
        return counts

    def nearbyKmers(self, seed, theta=-np.inf, window=6, k=3,
                    includePosition=False):
        s = self.sa.loc[range(*self.findKmer(seed))].copy()
        s = s.iloc[(self.scores.iloc[self.sourceBlock(s)] >= theta).values]
        if not "__contains__" in dir(window):
            window = [len(seed), window + len(seed)]
        strings = self.kmers(s, k=window[1], k0=window[0],
                             sanitize=False).loc[s]
        strLens = strings.str.len()
        if includePosition:
            out = pd.concat([str(i) + ' ' +
                             strings.loc[strLens >= i+k].str[i:(i+k)]
                             for i in range(strLens.max()-k+1)])
            out = out.value_counts()
            out = pd.DataFrame({
                "pos" : window[0] + out.index.str.replace(r' .*',
                                                          '').astype(int),
                "k" : k,
                "kmer" : out.index.str.replace(r'^\d+ ', ''),
                "count" : out
            }, index=out.index)
            return out
        else:
            return pd.concat([strings.loc[strLens >= i+k].str[i:(i+k)]
                              for i in range(strLens.max()-k+1)])

    def rankNearbyKmers(self, seed, theta=-np.inf, window=10,
                        ks=range(3, 7), minCount=2, bg=None, maxP=None):
        if bg is None:
            bg = pd.Series(np.ones(4), index=['A', 'C', 'G', 'T']) / 4.0
        if not "__contains__" in dir(window):
            window = [len(seed), window + len(seed)]
        if "__contains__" not in dir(ks):
            ks = [ks]
        ks = [k for k in ks if k <= window[1]-window[0]]
        kmers = pd.concat([self.nearbyKmers(seed, theta, window, k, True)
                           for k in ks])
        nSeqs = kmers.loc[(kmers['pos'] == kmers['pos'].min()) &
                          (kmers['k'] == kmers['k'].min()),
                          'count'].sum()
        kmers = kmers.loc[kmers['count'] >= minCount]
        kmers = kmers.loc[~kmers['kmer'].str.match(r".*\$")]
        kmers['logp'] = [
            stats.binom(
                n = nSeqs, 
                p = np.prod([bg[letter] for letter in kmers.ix[i, "kmer"]])
            ).logsf(kmers.ix[i, 'count'] - 1)
            for i in range(kmers.shape[0])
        ] / np.log(10)
        if kmers.shape[0] > 0 and maxP is not None:
            kmers = kmers.loc[kmers['logp'] <= np.log10(maxP)]
        return kmers.sort_values('logp')

    def clusterKmers(kmers, d=None):
        starcode = "starcode -s --print-clusters"
        if d is not None:
            starcode += " -d " + str(d)
        scOut = check_output(
            starcode,
            input = "\n".join([km + "\t" + str(len(kmers)-i)
                               for i, km in enumerate(kmers)]),
            universal_newlines = True,
            shell = True
        )
        clusters = {s.split("\t")[0] : s.split("\t")[2].split(",")
                    for s in scOut.split("\n") if s != ""}
        return clusters
        
    def alignCluster(self, kmers, theta=None,
                     start=-3, length=20,
                     minCount=None, minFrac=0.5,
                     msaFile=None, logoFile=None):
        seqFile = tempfile.NamedTemporaryFile()
        if msaFile is None:
            msaFile = tempfile.NamedTemporaryFile()
            msaFilename = msaFile.name
        else:
            msaFilename = msaFile
        kmerBounds = [self.findKmer(km) for km in kmers]
        clustEsses = [self.sa[kmb[0]:kmb[1]] for kmb in kmerBounds]
        clustEsses= pd.concat(clustEsses)
        if theta is not None:
            clustEsses = clustEsses.iloc[(self.scores.iloc[
                    self.sourceBlock(clustEsses)] >= theta).values]
        clustSeqs = self.kmers(clustEsses+start, k=length)
        clustSeqs = clustSeqs.loc[~clustSeqs.str.match(r".*\$")]
        with open(seqFile.name, "w") as outhandle:
            for seqName in clustSeqs.index:
                outhandle.write(">" + str(seqName) + "\n")
                outhandle.write(str(clustSeqs.loc[seqName]) + "\n")
        check_output("mafft " + seqFile.name + " > " + msaFilename, shell=True)
        seqFile.close()
        msaFasta = Fasta(msaFilename)
        seqNames = list(msaFasta.keys())
        alignMatrix = pd.DataFrame([list(str(msaFasta[seqName]).upper())
                                    for seqName in seqNames])
        alphabet = sorted(list(set(alignMatrix.values.flatten()) - set('-')))
        counts = pd.DataFrame(np.zeros((len(msaFasta[seqNames[0]]),
                                        len(alphabet))),
                              columns = alphabet)
        counts.index = range(counts.shape[0])
        for letter in alphabet:
            counts[letter] = (alignMatrix == letter).sum(axis=0)
        pwm = Sarks.pwmClean(counts, minCount, minFrac)
        if logoFile is not None:
            subFile = tempfile.NamedTemporaryFile()
            check_output("fasta-substring.py " + msaFilename +
                         " " + str(pwm.index.min()+1) +
                         " " + str(pwm.index.max()+1) +
                         " " + subFile.name,
                         shell=True)
            check_output("weblogo -A dna -D fasta -F pdf < "
                         + subFile.name + " > " + logoFile,
                         shell=True)
            subFile.close()
        return pwm

    def pwmClean(counts, minCount=None, minFrac=0.5):
        countSums = counts.sum(axis=1)
        if minCount is None:
            minCount = countSums.max() * 0.95
        pwm = counts.divide(countSums, axis=0)
        pwmMaxes = pwm.max(axis=1)
        start = 0
        while countSums.loc[start] < minCount or pwmMaxes.loc[start] < minFrac:
            start += 1
        end = counts.shape[0]
        while countSums.loc[end-1] < minCount or pwmMaxes.loc[end-1] < minFrac:
            end -= 1
        return pwm.loc[range(start, end)]

    def permute(self, perm):
        permutedScores = pd.Series(self.scores.iloc[perm].values,
                                   index = self.scores.index)
        permutedSarks = Sarks(self.inFasta, self.catFasta, permutedScores,
                              self.suffixArrayFile,
                              self.halfWindow, self.minScore,
                              self.maxScore, self.maxRun,
                              self.seqs, self.catSeq, self.bounds, self.sa,
                              self.windGini, self.spatialLength, self.spatGini)
        return permutedSarks

    def permutationDistribution(self, reps, k=12,
                                quantiles=None,
                                spatialLength=None,
                                minGini=None,
                                minSpatialGini=None,
                                seed=None):
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        for r in range(reps):
            rperm = np.random.choice(sl, sl, replace=False)
            rsarks = self.permute(rperm)
            rpos = rsarks.filter(minGini = minGini,
                                 spatialLength = spatialLength,
                                 minSpatialGini = minSpatialGini,
                                 k = k,
                                 prune = False,
                                 deduplicate = False)
            permResults.append({'windowed' : rsarks.windowed.loc[rpos]})
            if spatialLength is None and self.spatialLength is not None:
                spatialLength = self.spatialLength
            if spatialLength is not None:
                permResults[-1]['spatial_windowed'] = \
                        rsarks.spatialWindow(spatialLength).loc[rpos].dropna()
            if quantiles is not None:
                permResults[-1]['windowed'] = np.percentile(
                        permResults[-1]['windowed'], 100*quantiles)
                if spatialLength is not None:
                    permResults[-1]['spatial_windowed'] = np.percentile(
                            permResults[-1]['spatial_windowed'], 100*quantiles)
        return permResults

    def permutationTest(self, theta, reps,
                        k=12,
                        spatialLength=None, spatialTheta=None,
                        nearbyPars=None,
                        minGini=None, minSpatialGini=None,
                        seed=None):
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        for r in range(reps):
            rperm = np.random.choice(sl, sl, replace=False)
            rsarks = self.permute(rperm)
            rpos = rsarks.filter(minGini = minGini,
                                 theta = theta,
                                 spatialLength = spatialLength,
                                 minSpatialGini = minSpatialGini,
                                 spatialTheta = spatialTheta,
                                 nearbyPars = nearbyPars,
                                 k = k,
                                 prune = False,
                                 deduplicate = False)
            if len(rpos) == 0:
                permResults.append(pd.DataFrame({
                    'i':[], 's':[], 'kmer':[], 'khat':[],
                    'block':[], 'wi':[], 'gini':[],
                    'score':[], 'windowed':[]
                })[['i', 's', 'kmer', 'khat', 'block', 'wi',
                    'gini', 'score', 'windowed']])
                if spatialTheta is not None:
                    permResults[-1]['spatial_windowed'] = []
            else:
                rout = rsarks.subtable(rpos, k=k)
                rout = rout.loc[~rout['kmer'].str.match(r".*\$")]
                permResults.append(rout)
        return permResults
