from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import editdistance
import intervaltree
from io import StringIO
import numpy as np
import os
import pandas as pd
from pyfaidx import Fasta
import re
from scipy import stats
from subprocess import check_output
import tempfile


class Sarks(object):
    """
    Sarks class implements suffix array kernel smoothing for de novo correlative motif discovery, conducted in several steps:

    1. initialize Sarks object (requires fasta file and pandas Series containing scores, as well as specification of smoothing window size),
    2. assessment of distribution of scores under null-hypothesis using permutationDistribution specifying Gini impurity cutoffs (minGini and/or minSpatialGini),
    3. call to peaks method specifying filter parameters (e.g., minGini and smoothed score cutoff theta),
    4. extraction of kmer table for top peaks using subtable method (followed by optional use of extendKmers if desired),
    5. clustering of extracted kmers using clusterKmers."""

    def __init__(self, inFasta, catFasta, scores,
                 suffixArrayFile=None,
                 halfWindow=250,
                 seqs=None, catSeq=None, bounds=None, sa=None,
                 windGini=None, spatialLength=None, spatGini=None,
                 regenerateSuffixArray=False):
        """
        Construction of Sarks object requires, at minimum:

        :param inFasta: specification of fasta file containing sequences to be analyzed
        :param catFasta: file name to be used for fasta file containing concatenated sequence produced by Sarks
        :param scores: pandas Series object with index matching the sequence ids in inFasta
        :param halfWindow: half-width of smoothing window
        """
        self.inFasta = inFasta
        self.catFasta = catFasta
        self.suffixArrayFile = suffixArrayFile
        self.halfWindow = halfWindow
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
        """
        Concatenate seqs together and write to catFasta
        
        :param catFasta: name of file to write concatenated sequence to
        :param regenerateSuffixArray: if True, overwrite any existing concatenated sequence file
        """
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
        self.catSeq = Fasta(catFasta)['concatenated_promoters']
        self.bounds = np.cumsum([len(self.seqs[t]) + 1
                                 for t in self.transcripts])

    def calcSuffixArray(self, catFasta, suffixArrayFile,
                        regenerateSuffixArray=False):
        """
        Construct suffix array for concatenated sequence held in catFasta

        :param catFasta: name of file containing concatenated sequence
        :param suffixArrayFile: name of file to write suffix array data to (or load from if already exists and regenerateSuffixArray is False)
        :param regenerateSuffixArray: if True, overwrite any existing suffix array file
        """
        if suffixArrayFile is None:
            suffixArrayFile = re.sub(r"\..*", "", catFasta) + "_sa.txt"
        # if regenerateSuffixArray or not os.path.exists(suffixArrayFile + ".gz"):
        #     check_output(
        #         "suffix-array " + catFasta +
        #         " | grep -v '>' > " + suffixArrayFile,
        #         shell = True
        #     )
        #     check_output("gzip -f " + suffixArrayFile, shell=True)
        # self.sa = pd.Series(np.loadtxt(suffixArrayFile + ".gz").astype(int))
        self.sa = pd.read_csv(StringIO(check_output(
            'suffix-array ' + catFasta + ' | grep -v ">"',
            shell = True
        ).decode('utf-8')), index_col=None, header=None).iloc[:, 0].astype(int)
        ## need to remove extra end-of-string position
        self.sa = self.sa.loc[self.sa < len(self.catSeq)]
        
    def sourceBlock(self, s):
        """
        Map suffix array *value* s back to numeric index of source block

        :param s: suffix array value (position of suffix within concatenated sequence)
        :returns: numeric index of source block in the Series scores
        """
        return np.searchsorted(self.bounds, s+1)

    def sourceScore(self, s):
        """
        Map suffix array *value* s back to score of block from which s is derived

        :param s: suffix array value (position of suffix within concatenated sequence)
        :returns: score of block from which s is derived
        """
        srcBlk = self.sourceBlock(s)
        out = pd.DataFrame({
            'block' : srcBlk,
            # 'score' : self.scores.iloc[srcBlk].values
            'score' : self.scores.values[srcBlk]
        }, index=s)
        out['block'] = out['block'].astype(int)
        return out[['block', 'score']]

    def windowGini(self, recalculate=False):
        """
        Calculate Gini impurities for all positions
        
        :param recalculate: if False and windGini already set simply return existing values
        :returns: pandas Series indexed by suffix array *value* s (position within concatenated sequence); also stored internally as windGini
        """
        if "windGini" not in dir(self) or recalculate:
            block = pd.Series(self.sourceBlock(self.sa), index=self.sa)
            btmp = tempfile.NamedTemporaryFile()
            np.savetxt(btmp.name, block.values, '%d')
            gtmp = tempfile.NamedTemporaryFile()
            check_output(
                "windginiimp " + str(self.halfWindow) + ' ' + btmp.name + \
                ' > ' + gtmp.name,
                shell = True
            )
            self.windGini = pd.Series(np.loadtxt(gtmp.name))
            btmp.close()
            gtmp.close()
            # saWindowCenters = range(
            #     int(self.halfWindow),
            #     int(self.halfWindow) + len(self.windGini)
            # )
            # self.windGini.index = self.sa.iloc[saWindowCenters]
            self.windGini.index = self.sa.values[
                    int(self.halfWindow):(int(self.halfWindow)+len(self.windGini))]
        return self.windGini

    def spatialGini(self, recalculate=False):
        """
        Calculate Gini impurities, averaged over a spatial window of length spatialLength, for all positions

        :param recalculate: if False and windGini already set simply return existing values
        :returns: pandas Series indexed by suffix array *value* s (position within concatenated sequence) of start of spatial window; also stored internally as spatGini
        """
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
                    sortedGini.index.values[0:len(self.spatGini)]
                    # sortedGini.iloc[range(len(self.spatGini))].index
        return self.spatGini
    
    def window(self, halfWindow=None, recalculate=False):
        """
        Access to sorted-suffix-smoothed scores yhat

        :param halfWindow: number of suffixes before and after central suffix in window to smooth over
        :param recalculate: if False and windGini already set simply return existing values
        :returns: pandas Series indexed by suffix array *value* s (position within concatenated sequence)
        """
        if halfWindow is None:
            halfWindow = self.halfWindow
        if halfWindow != self.halfWindow or recalculate:
            if halfWindow != self.halfWindow:
                self.halfWindow = halfWindow
                self.windowGini(recalculate=True)
                self.spatGini = None
                if 'spatialWindowed' in dir(self):
                    del self.spatialWindowed
            saScores = self.sourceScore(self.sa)['score']
            saCumul = saScores.cumsum()
            windowSize = (2 * halfWindow) + 1
            self.windowed = saCumul.iloc[(windowSize-1):] - \
                            np.concatenate([
                                np.zeros(1),
                                saCumul.iloc[0:(-windowSize)].values
                            ])
            self.windowed /= windowSize
            # saWindowCenters = np.arange(
            #     int(halfWindow),
            #     int(halfWindow) + len(self.windowed)
            # )
            # self.windowed.index = saScores.iloc[saWindowCenters].index
            self.windowed.index = saScores.index.values[
                    int(halfWindow):(int(halfWindow)+len(self.windowed))]
        return self.windowed

    def spatialWindow(self, length):
        """
        Access to spatially-smoothed scores yhathat
        
        :param length: spatial distance to smooth over
        :returns: pandas Series indexed by suffix array *value* s (position within concatenated sequence)
        """
        if 'spatialLength' not in dir(self) or \
           int(length) != self.spatialLength or \
           'spatialWindowed' not in dir(self):
            self.spatGini = None
            self.spatialLength = int(length)
            windowed = self.window().copy()
            windowed.sort_index(inplace=True)
            windCumul = windowed.cumsum()
            self.spatialWindowed = windCumul.iloc[(int(length)-1):] - \
                                   np.concatenate([
                                       np.zeros(1),
                                       windCumul.iloc[0:(-int(length))].values
                                   ])
            self.spatialWindowed /= int(length)
            self.spatialWindowed.index = \
                    windowed.index.values[0:len(self.spatialWindowed)]
                    # windowed.iloc[range(len(self.spatialWindowed))].index
        return self.spatialWindowed

    def i2s(self, i):
        """
        Get suffix array *value* s corresponding to given suffix array *index* i

        :param i: suffix array index (position of suffix in lexicographically sorted list of suffixes) (int)
        :returns: suffix array value s (position with concatenated sequence) (int)
        """
        return int(self.sa.loc[i])

    def s2i(self, s):
        """
        Get suffix array *index* i corresponding to given suffix array *value* s

        :param s: suffix array value (position with concatenated sequence) (int)
        :returns: suffix array index (position of suffix in lexicographically sorted list of suffixes) (int)
        """
        return int(self.sa.loc[self.sa == s].index[0])

    def findKmer(self, kmer):
        """
        Find interval [imin, imax) such that all suffix array index positions imin <= i < imax correspond to suffixes beginning with string kmer

        :param kmer: string representing the common prefix to all suffixes of concatenated sequence stored in catFasta whose suffix array index is >= imin and < imax
        :returns: tuple (imin, imax)
        """
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
        """
        Use intervaltree to remove any nested intervals

        :param intervals: list of tuples representing intervals to prune
        :returns: list of tuple representations of intervals remaining after pruning
        """
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
              deduplicate=False):
        """
        Identify local score maxima when sorted by suffix array index i

        :param theta: minimum smoothed score value for peak to be reported (float)
        :param prune: whether to apply pruneIntervals to resulting kmer intervals (bool)
        :param spatialLength: length of spatial window to apply (if any) (int)
        :param spatialTheta: minimum spatially smoothed score value for peak to be reported (if any; requires spatialLength to be specified) (float)
        :param k: maximum kmer length to consider (int)
        :param minK: minimum kmer length to report (int)
        :param minGini: minimum Gini impurity value for suffix position to be reported as peak (float)
        :param minSpatialGini: minimum spatially-averaged Gini impurity value for suffix position to be reported as peak (float)
        :param deduplicate: if True, report only highest scoring position corresponding to any given kmer (bool)
        :returns: suffix array *values* s corresponding to desired peak positions (pandas Series)
        """
        scoreDiffB = self.windowed.sort_index().diff()
        scoreDiffF = scoreDiffB.copy().iloc[1:].copy()
        scoreDiffF.index = scoreDiffB.index[:-1]
        scoreDiffB = scoreDiffB.iloc[:-1]
        pos = self.filter(minGini=minGini, minSpatialGini=minSpatialGini,
                          minK=minK,
                          theta=theta,
                          spatialLength=spatialLength, spatialTheta=spatialTheta,
                          nearbyPars=None,
                          k=k, prune=prune, deduplicate=deduplicate)
        pos = pos.loc[(scoreDiffB.loc[pos] >= 0) &
                      (scoreDiffF.loc[pos] <= 0)]
        return self.sa.loc[self.sa.isin(set(pos))]

    def filter(self, minGini=None, minSpatialGini=None, minK=None,
               theta=None, spatialLength=None, spatialTheta=None,
               nearbyPars=None, k=12, prune=True, deduplicate=True):
        """
        Return suffix array *values* s satisfying specified filters
        
        :param minGini: minimum Gini impurity value for suffix position to be reported as peak (float)
        :param minSpatialGini: minimum spatially-averaged Gini impurity value for suffix position to be reported as peak (float)
        :param minK: minimum kmer length to report (int)
        :param theta: minimum smoothed score value for peak to be reported (float)
        :param spatialLength: length of spatial window to apply (if any) (int)
        :param spatialTheta: minimum spatially smoothed score value for peak to be reported (if any; requires spatialLength to be specified) (float)
        :param nearbyPars: args for nearbyKmers (dict)
        :param k: maximum kmer length to consider (int)
        :param prune: whether to apply pruneIntervals to resulting kmer intervals (bool)
        :param deduplicate: if True, report only highest scoring position corresponding to any given kmer (bool)
        :returns: suffix array *values* s satisfying specified filters (pandas Series)
        """
        pos = pd.Series(self.windowed.index)
        pos.index = pos
        posSubtable = None
        if minGini is not None:
            if minGini >= 1:
                minGini = 1 - (1 - self.windGini.median()) * minGini
            pos = pos.loc[self.windGini[pos] >= minGini]
        if theta is not None:
            pos = pos.loc[(self.windowed.loc[pos] >= theta)]
        if spatialLength is not None and spatialLength > 0:
            spat = self.spatialWindow(spatialLength)
            if minSpatialGini is not None:
                if minSpatialGini >= 1:
                    minSpatialGini = 1 - (1 - self.windGini.median()) * minSpatialGini
                pos = pos.loc[self.spatialGini()[pos] >= minSpatialGini]
            if spatialTheta is not None:
                pos = pos.loc[spat.loc[pos] >= spatialTheta]
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
        """
        Extend kmers when adding neighboring characters would result in another reported kmer string

        :param subtable: result of calling subtable method on filtered peaks (pandas DataFrame)
        :returns: subtable with extended kmer values in kmer column (pandas DataFrame)
        """
        subtable = subtable.copy()
        maxLen = subtable['kmer'].str.len().max()
        oldS = pd.Series(-1, index=subtable.index)
        while not (oldS == subtable['s']).all():
            oldS = subtable['s'].copy()
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

    def kmers(self, s, k=12, k0=0, sanitize=True):
        """
        Report (k-k0)-mers starting at positions s+k0 ending at s+k

        :param s: list/array/Series of suffix array *values* at which kmers are located
        :param k: int offset specifying how many characters downstream of s kmers should end
        :param k0: int offset specifying how many characters downstream of s kmers should start
        :param sanitize: if True, remove any kmers not of length (k-k0) (from end of concatenated sequence)
        :returns: (k-k0)-mers starting at positions s+k0 ending at s+k (pandas Series)
        """
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

    # @staticmethod
    # def prefixAgreeSum(kmer, esses):
    #     """
    #     Calculate average length of prefix agreement between kmer and suffixes specified by esses
        
    #     :param kmer: string to assess average length of maximum prefix agreement with
    #     :param esses: suffix array *values* (positions within concatenated sequence) of suffixes to compare with prefix
    #     :returns: average length of prefix agreement between kmer and suffixes specified by esses (float)
    #     """
    #     esses = pd.Series(esses)
    #     out = 0
    #     aliveIndices = esses.index
    #     for i, letter in enumerate(str(kmer)):
    #         iLetters = esses[aliveIndices].str[i].dropna()
    #         aliveIndices = esses[iLetters.index][iLetters == letter].index
    #         out += len(aliveIndices)
    #     return out / len(esses)

    def prefixAgreeSum(self, i, halfWindow=None, kmax=12):
        """
        Calculate average length of prefix agreement between kmer specified by suffix array *index* i and suffixes in smoothing window

        :param i: suffix array *index* at center of smoothing window (int)
        :param halfWindow: half-width of smoothing window (int)
        :param kmax: maximum kmer length to consider (int)
        """
        if halfWindow is None:
            halfWindow = self.halfWindow
        kmer = self.kmers([self.sa[i]]).iloc[0]
        agreeSum = 0
        wstart, wend = i-halfWindow, i+halfWindow+1
        kstart, kend = i, i+1
        for k in range(kmax, 0, -1):
            kwin = list(self.findKmer(kmer[0:k]))
            kwin[0] = max(wstart, kwin[0])
            kwin[1] = min(wend, kwin[1])
            agreeSum += (k * ((kstart - kwin[0]) + (kwin[1] - kend)))
            kstart, kend = kwin
        return agreeSum / (2 * halfWindow)

    def subtable(self, s, k=12):
        """
        Suffix array table for suffix array *values* s

        :param s: pandas Series containing suffix array values (positions within concatenated sequence)
        :param k: maximum kmer length to consider (int)
        :returns: pandas DataFrame containing info, including conserved kmers, for positions specified by s
        """
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
            "khat" : [self.prefixAgreeSum(i, self.halfWindow, k)
                      for i in eyes],
            # "khat" : [Sarks.prefixAgreeSum(
            #     km,
            #     self.kmers(self.sa[
            #         list(range(i-self.halfWindow, i)) +
            #         list(range(i+1, i+self.halfWindow+1))
            #     ], k, 0, False)) for i, km in zip(eyes, kmers)],
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
        """
        Suffix array table for suffix array *values* s including kmers of pre-specified length k

        :param s: pandas Series containing suffix array values (positions within concatenated sequence)
        :param k: kmer length (int)
        :returns: pandas DataFrame containing info, including constant length kmers, for positions specified by s
        """
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

    def nearbyKmers(self, seed, theta=-np.inf, window=6, k=3,
                    includePosition=False):
        """
        Count kmers downstream of specified seed sequence

        :param seed: string representing head kmer for which enriched downstream tail kmers are sought
        :param theta: minimum smoothed score for suffixes beginning with seed to be included in analysis (float)
        :param window: maximum distance downstream of (end of) seed to consider (int)
        :param k: minimum length of enriched kmer to analyze (int)
        :param includePosition: format output as pandas DataFrame (otherwise returns Series)
        :returns: either pandas Series or DataFrame counting kmers downstream of specified seed sequence
        """
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
        """
        Identify potentially enriched kmers downstream of specified seed sequence

        :param seed: string representing head kmer for which enriched downstream tail kmers are sought
        :param theta: minimum smoothed score for suffixes beginning with seed to be included in analysis (float)
        :param window: maximum distance downstream of (end of) seed to consider (int)
        :param ks: tuple indicating (minimum, maximum+1) length of kmers to consider (tuple of ints)
        :param minCount: minimum number of occurrences of kmer for consideration (int)
        :param bg: null model probability of each possible character (pandas Series indexed by characters in sequence alphabet)
        :param maxP: maximum p-value allowed to retain kmers in output (float)
        :returns: table of info for potentially enriched downstream kmers (pandas DataFrame)
        """
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

    @staticmethod
    def clusterKmers(kmers, d=None):
        """
        Uses starcode (-s --print-clusters -d d) to cluster specified kmers

        :param kmers: list/array/Series of kmers to cluster (higher priority kmers first)
        :param d: maximum edit distance between kmer and cluster center for inclusion (int)
        :returns: dict mapping cluster center to list of kmers included in cluster
        """
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

    @staticmethod
    def revCompFilter(kmers, d=None):
        kmers = list(pd.Series(list(kmers)).str.replace('\$.*$', ''))
        def revComp(s):
            comp = {'a':'t', 'A':'T', 'c':'g', 'C':'G', 'n':'n', 'N':'N',
                    'g':'c', 'G':'C', 't':'a', 'T':'A'}
            return ''.join([comp[s[i]] for i in range(len(s)-1, -1, -1)])
        retain = set()
        for i, km in enumerate(kmers):
            rkm = revComp(km)
            for j in range(i+1, len(kmers)):
                if editdistance.eval(rkm, kmers[j]) <= d:
                    retain.add(km)
                    retain.add(kmers[j])
        return sorted(list(retain))

    # def estimateSignalToNoise(self, halfWindows, minGinis):
    #     out = []
    #     scoreMean = self.scores.mean()
    #     scoreVar = self.scores.var(ddof=0)
    #     for halfWindow in halfWindows:
    #         windowed = self.window(halfWindow)
    #         for minGini in minGinis:
    #             acceptGini = (self.windGini >= minGini)
    #             sig = windowed.loc[acceptGini] - scoreMean
    #             estPermVar = scoreVar *\
    #                          (1 - self.windGini.loc[acceptGini])
    #             out.append((halfWindow,
    #                         minGini,
    #                         estPermVar.max(),
    #                         sig.max() / np.sqrt(estPermVar.max())))
    #     return pd.DataFrame(out, columns=[
    #             'half_window', 'min_gini', 'est_perm_var', 'signal_to_noise'])
        
    def permute(self, perm):
        """
        Generate new Sarks object with scores permuted by perm

        :param perm: array representing score permutation (e.g., output from np.random.choice)
        :returns: new Sarks object with scores permuted by perm; where possible, shares member variales with self
        """
        permutedScores = pd.Series(self.scores.iloc[perm].values,
                                   index = self.scores.index)
        permutedSarks = Sarks(inFasta = self.inFasta,
                              catFasta = self.catFasta,
                              scores = permutedScores,
                              suffixArrayFile = self.suffixArrayFile,
                              halfWindow = self.halfWindow,
                              seqs = self.seqs,
                              catSeq = self.catSeq,
                              bounds = self.bounds,
                              sa = self.sa,
                              windGini = self.windGini,
                              spatialLength = self.spatialLength,
                              spatGini = self.spatGini,
                              regenerateSuffixArray = False)
        return permutedSarks

    def permutationDistribution(self, reps, k=12,
                                quantiles=np.array([0, 0.5, 0.99, 0.9999,
                                                    0.999999, 1]),
                                spatialLength=None,
                                minGini=None,
                                minSpatialGini=None,
                                seed=None,
                                permutations=None):
        """
        Estimate null distribution of smoothed (and spatially-smoothed, if desired) scores passing specified filters

        :param reps: how many repetitions of random permutation to do (int)
        :param k: maximum kmer length to consider (int)
        :param quantiles: quantiles (between 0 and 1) of distribution to report (array or Series)
        :param spatialLength: length of spatial window for spatial smoothing (int)
        :param minGini: minimum Gini impurity value for suffix position to be included (float)
        :param minSpatialGini: minimum spatially-averaged Gini impurity value for suffix position to be included (float)
        :param seed: seed for random number generator (int)
        :param permutations: reps-in-rows, entries-in-columns table of permutations to apply (pandas DataFrame)
        :returns: returns list of dicts (keyed by 'windowed' and, if applicable, 'spatial_windowed') containing requested quantiles of filtered smoothed scores
        """
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        for r in range(reps):
            if permutations is not None:
                rperm = permutations.values[r, :]
            else:
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
            if spatialLength is not None and spatialLength > 0:
                permResults[-1]['spatial_windowed'] = \
                        rsarks.spatialWindow(spatialLength).loc[rpos].dropna()
            if quantiles is not None:
                permResults[-1]['windowed'] = np.percentile(
                        permResults[-1]['windowed'], list(100*quantiles))
                if spatialLength is not None and spatialLength > 0:
                    permResults[-1]['spatial_windowed'] = np.percentile(
                            permResults[-1]['spatial_windowed'], list(100*quantiles))
        return permResults

    def permutationDistributionMultiFilter(self, reps, k=12,
                                           quantiles=np.array([0, 0.5, 0.99,
                                                               0.9999,
                                                               0.999999, 1]),
                                           filters = None,
                                           seed = None,
                                           permutations = None):
        """
        Estimate null distribution of smoothed (and spatially-smoothed, if desired) scores passing specified filters

        :param reps: how many repetitions of random permutation to do (int)
        :param k: maximum kmer length to consider (int)
        :param quantiles: quantiles (between 0 and 1) of distribution to report (array or Series)
        :param filters: table providing sets of spatialLength, minGini, and minSpatialGini values to assess (pandas DataFrame)
        :param seed: seed for random number generator (int)
        :param permutations: reps-in-rows, entries-in-columns table of permutations to apply (pandas DataFrame)
        :returns: returns dict (keyed by 'windowed' and 'spatial_windowed') containing pandas DataFrames in turn containing requested quantiles of filtered smoothed scores
        """
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        out = {'windowed' : [], 'spatial_windowed' : []}        
        for r in range(reps):
            if permutations is not None:
                rperm = permutations.values[r, :]
            else:
                rperm = np.random.choice(sl, sl, replace=False)
            rsarks = self.permute(rperm)
            filts = filters.sort_values('spatialLength')
            for filtIdx in filts.index:
                filt = filts.loc[filtIdx]
                rpos = rsarks.filter(minGini = filt['minGini'],
                                     spatialLength = filt['spatialLength'],
                                     minSpatialGini = filt['minSpatialGini'],
                                     k = k,
                                     prune = False,
                                     deduplicate = False)
                windQuants = np.percentile(rsarks.windowed.loc[rpos],
                                           list(100 * quantiles))
                windQuants = pd.DataFrame(
                    windQuants,
                    index = quantiles,
                    columns = [len(out['windowed'])]
                ).transpose()
                windQuants['rep'] = r
                windQuants['halfWindow'] = self.halfWindow
                windQuants['spatialLength'] = filt['spatialLength']
                windQuants['minGini'] = filt['minGini']
                windQuants['minSpatialGini'] = filt['minSpatialGini']
                out['windowed'].append(windQuants)
                if filt['spatialLength'] > 0:
                    spatQuants = np.percentile(
                        rsarks.spatialWindow(filt['spatialLength'])\
                        .loc[rpos].dropna(),
                        list(100 * quantiles)
                    )
                    spatQuants = pd.DataFrame(
                        spatQuants,
                        index = quantiles,
                        columns = [len(out['spatial_windowed'])]
                    ).transpose()
                    spatQuants['rep'] = r
                    spatQuants['halfWindow'] = self.halfWindow
                    spatQuants['spatialLength'] = filt['spatialLength']
                    spatQuants['minGini'] = filt['minGini']
                    spatQuants['minSpatialGini'] = filt['minSpatialGini']
                    out['spatial_windowed'].append(spatQuants)
        out['windowed'] = pd.concat(out['windowed'])
        if len(out['spatial_windowed']) > 0:
            out['spatial_windowed'] = pd.concat(out['spatial_windowed'])
        else:
            out['spatial_windowed'] = None
        return out

    def permutationTestMultiFilter(self, reps, k=12,
                                   filters = None,
                                   seed = None,
                                   permutations = None):
        """
        Estimate null distribution of smoothed (and spatially-smoothed, if desired) scores passing specified filters

        :param reps: how many repetitions of random permutation to do (int)
        :param k: maximum kmer length to consider (int)
        :param filters: table providing sets of spatialLength, minGini, minSpatialGini, theta, and spatialTheta values to assess (pandas DataFrame)
        :param seed: seed for random number generator (int)
        :param permutations: reps-in-rows, entries-in-columns table of permutations to apply (pandas DataFrame)
        :returns: returns pandas DataFrame indicating count of positive hits for each combination of parameters specified by a single row of filters
        """
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        out = []
        for r in range(reps):
            if permutations is not None:
                rperm = permutations.values[r, :]
            else:
                rperm = np.random.choice(sl, sl, replace=False)
            rsarks = self.permute(rperm)
            filts = filters.sort_values(['halfWindow', 'spatialLength'])
            for filtIdx in filts.index:
                filt = filts.loc[filtIdx]
                rsarks.window(int(filt['halfWindow']))
                spatTheta = filt['spatialTheta']
                if np.isnan(spatTheta):
                    spatTheta = None
                rpos = rsarks.filter(minGini = filt['minGini'],
                                     spatialLength = filt['spatialLength'],
                                     minSpatialGini = filt['minSpatialGini'],
                                     theta = filt['theta'],
                                     spatialTheta = spatTheta,
                                     k = filt['kmax'] if 'kmax' in filt.index else k,
                                     prune = False,
                                     deduplicate = False)
                if 'd' in filt.index:
                    rkm = rsarks.kmers(rpos,
                                       k = filt['k'] if 'k' in filt.index else k)
                    rpos = rpos.loc[rkm.isin(revCompFilter(rkm, d=filt['d']))]
                rout = filt.copy()
                rout['rep'] = r
                rout['count'] = len(rpos)
                out.append(rout)
        return pd.concat(out, axis=1, ignore_index=True).transpose()
    
    def permutationTest(self, theta, reps,
                        k=12,
                        spatialLength=None, spatialTheta=None,
                        nearbyPars=None,
                        minGini=None, minSpatialGini=None,
                        seed=None, permutations=None):
        """
        Return filtered peaks, if any, obtained under permutation of input scores

        :param theta: minimum smoothed score value for peak to be reported (float)
        :param reps: how many repetitions of random permutation to do (int)
        :param k: maximum kmer length to consider (int)
        :param spatialLength: length of spatial window for spatial smoothing (int)
        :param spatialTheta: minimum spatially smoothed score value for peak to be reported (if any; requires spatialLength to be specified) (float)
        :param nearbyPars: args for nearbyKmers (dict)
        :param minGini: minimum Gini impurity value for suffix position to be included (float)
        :param minSpatialGini: minimum spatially-averaged Gini impurity value for suffix position to be included (float)
        :param seed: seed for random number generator (int)
        :param permutations: reps-in-rows, entries-in-columns table of permutations to apply (pandas DataFrame)
        :returns: returns suffix array subtables corresponding to peaks detected under each randonmly-generated permutation (list)
        """
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        for r in range(reps):
            if permutations is not None:
                rperm = permutations.values[r, :]
            else:
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

    def multiWindowPermute(
            self,
            halfWindows, filters, reps, k=12,
            quantiles=np.array([0, 0.5, 0.99, 0.9999, 0.999999, 1]),
            nsigma=None,
            seed=None, permutations=None, assess='distribution'):
        if seed is not None:
            np.random.seed(seed)
        permResults = []
        sl = len(self.scores)
        if permutations is None:
            permutations = pd.DataFrame([
                np.random.choice(sl, sl, replace=False)
                for r in range(reps)
            ])
        out = []
        for halfWindow in halfWindows:
            self.window(int(halfWindow))
            if assess == 'distribution':
                out.append(self.permutationDistributionMultiFilter(
                    reps=reps, k=k, quantiles=quantiles,
                    filters=filters, seed=seed, permutations=permutations
                ))
            elif assess == 'test':
                out.append(self.permutationTestMultiFilter(
                    reps=reps, k=k,
                    filters=filters, seed=seed, permutations=permutations
                ))
        if not isinstance(out[0], pd.DataFrame):
            preout = out
            out = {'windowed' : pd.concat([x['windowed'] for x in preout])}
            outSpat = [x['spatial_windowed']
                       for x in preout
                       if 'spatial_windowed' in x
                       and x['spatial_windowed'] is not None]
            if len(outSpat) > 0:
                out['spatial_windowed'] = pd.concat(outSpat)
            else:
                out['spatial_windowed'] = out['windowed'].copy().iloc[[]]
            if nsigma is not None:
                groupCols = ['halfWindow', 'spatialLength',
                             'minGini', 'minSpatialGini']
                theta = out['windowed'][[1.0] + groupCols]\
                        .groupby(groupCols)\
                        .agg(lambda x: x.mean() + nsigma*x.std(ddof=1))
                theta.columns = ['theta']
                spatTheta = out['spatial_windowed'][[1.0] + groupCols]\
                            .groupby(groupCols)\
                            .agg(lambda x: x.mean() + nsigma*x.std(ddof=1))
                if spatTheta.shape[0] > 0:
                    spatTheta.columns = ['spatialTheta']
                    theta['spatialTheta'] = spatTheta.loc[theta.index,
                                                          'spatialTheta']
                    theta.loc[~theta['spatialTheta'].isnull(), 'theta'] = -np.inf
                else:
                    theta['spatialTheta'] = -np.inf
                theta['kmax'] = k
                theta = theta[['kmax', 'theta', 'spatialTheta']]
                out['theta'] = theta.reset_index()
                # out['theta']['spatialLength'] =\
                #         out['theta']['spatialLength'].astype(int)
        else:
            out = pd.concat(out)
        return out

    def multiWindowPeaks(self, filters, k=12,
                         prune=False, extend=False, deduplicate=False):
        peakTables = OrderedDict()
        for filtIdx in filters.index:
            filt = filters.loc[filtIdx].copy()
            filtKey = tuple([(k, filt[k]) for k in filt.index])
            self.window(int(filt['halfWindow']))
            spatTheta = filt['spatialTheta']
            if np.isnan(spatTheta):
                spatTheta = None
            peaks = self.peaks(minGini = filt['minGini'],
                               spatialLength = int(filt['spatialLength']),
                               minSpatialGini = filt['minSpatialGini'],
                               theta = filt['theta'],
                               spatialTheta = spatTheta,
                               prune = prune,
                               deduplicate = deduplicate,
                               k = k)
            peakTable = self.subtable(peaks)
            if extend and peakTable.shape[0] > 0:
                peakTable = self.extendKmers(peakTable)
            peakTables[filtKey] = peakTable
        return peakTables
