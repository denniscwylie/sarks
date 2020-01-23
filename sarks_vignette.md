Suffix Array Kernel Smoothing for de novo discovery of correlative sequence motifs and multi-motif domains: the sarks package
------------------------------------------------------------------------------------------------------------------------------

This vignette describes the R interface sarks to the java implementation of the Suffix
Array Kernel Smoothing, or SArKS, algorithm for identification of correlative motifs and
multi-motif domains (MMDs).

If you use sarks in your work, please cite Wylie et al. [2019] (see [references](#references)).

Table of Contents
-----------------

   * [1 How to use sarks](#1-how-to-use-sarks)
      * [1.1 Starting sarks](#11-starting-sarks)
      * [1.2 Input](#12-input)
         * [1.2.1 SArKS run parameters](#121-sarks-run-parameters)
      * [1.3 Output](#13-output)
         * [1.3.1 k-mers](#131-k-mers)
         * [1.3.2 Reducing redundancy in peaks](#132-reducing-redundancy-in-peaks)
         * [1.3.3 Suffix arrays](#133-suffix-arrays)
         * [1.3.4 Window averaging (kernel smoothing)](#134-window-averaging-kernel-smoothing)
      * [1.4 Permutations and thresholds](#14-permutations-and-thresholds)
         * [1.4.1 Gini impurity filter](#141-gini-impurity-filter)
         * [1.4.2 Multithreading sarks permutational analyses](#142-multithreading-sarks-permutational-analyses)
      * [1.5 Spatial smoothing](#15-spatial-smoothing)
      * [1.6 Varying SArKS parameters](#16-varying-sarks-parameters)
      * [1.7 Clustering similar k-mers into broader motifs](#17-clustering-similar-k-mers-into-broader-motifs)
   * [2 How sarks works](#2-how-sarks-works) *omitted; see vignette pdf*
   * [3 Notation glossary](#3-notation-glossary) *omitted; see vignette pdf*
   * [References](#references)

# 1 How to use sarks

## 1.1 Starting sarks

Because sarks is implemented using rJava, and because the default setting for the java
virtual machine (JVM) heap space with rJava is quite low, you should initialize java with
increased heap size before loading sarks (or any other rJava-dependent R packages):

```R
options(java.parameters = '-Xmx8G') ## sets to 8 gigabytes: modify as needed
library(rJava)
.jinit()
## ---------------------------------------------------------------------
## once you've set the java.parameters and then called .jinit from rJava
## (making sure not to load any other rJava-dependent packages prior to
## initializing the JVM with .jinit!),
## you are now ready to load the sarks library:
## ----------------------------------------------------------------------
library(sarks)
```

## 1.2 Input

Aside from the specification of a few analysis parameters to be discussed below, sarks
requires two pieces of input data:
- sequences which may be specified as either a:
      FASTA formatted text file (may be gzipped) or, equivalently, a
      named character vector
- scores which may be specified as either a:
      named numeric vector (with names matching those of input sequences) or a
      tsv two-column tab-delimited file (containing header line), with columns
              1. containing names matching those of input sequences and
              2. containing numeric scores assigned to those sequences

For the purposes of this vignette, we’ll take the sequences to be the named character vector
simulatedSeqs included in the sarks package. simulatedSeqs is a character vector of
length 30—representing 30 different (fake) DNA sequences—each 250 characters in length:

Take a look at the first few characters of each sequence:

```R
vapply(simulatedSeqs, function(s) {paste0(substr(s, 1, 10), '...')}, '')
```
*out:*
```
              0                  1                   2                    3                 4
"GCGGAGGCTG..."    "CGTTGAATGT..."     "AGTCAGTTCT..."      "AGAGCTTCAG..."   "GTTTCTGCCC..."
              5                  6                   7                    8                 9
"CTAAGGGCGA..."    "ATTAGGTAAA..."     "GCTCGGAGGA..."      "TTCCTGCCTA..."   "ACAATCTGCG..."
             10                 11                  12                   13                14
"CCACAGCGTT..."    "TGACGACGCG..."     "GCGCACTAGC..."      "TCAAAGTAGG..."   "GGTACAATCA..."
             15                 16                  17                   18                19
"TATGACACCG..."    "CACTCGTATG..."     "TGGTCTCGAC..."      "GTCTCCCCGA..."   "TACGAGGCTC..."
             20                 21                  22                   23                24
"CGGACGCGTG..."    "GATGTGCCAT..."     "TGAAAGGAGA..."      "TAATGTAATG..."   "CATCGAGATG..."
             25                 26                  27                   28                29
"CATACTGAGA..."    "ACCAACAGTC..."     "GCACGACAAA..."      "GAAACAGAGG..."   "GTTGATCTCA..."
```

Notice that the names of the simulated sequences are just the (string representations of) the
numbers "0" through "29".

Now let’s look at the simulatedScores:

```R
simulatedScores
```
*out:*
```
 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1
26 27 28 29
 1  1  1  1
```

Let’s check that the names of the scores match those of the sequences:

```R
all(names(simulatedScores) == names(simulatedSeqs))
```
*out:*
```
[1] TRUE
```

Looking back at the scores themselves, note that for this simple illustrative example the
scores are just defined to be 0 for the first 20 sequences ("0"–"19") and 1 for the last 10
sequences ("20"–"29").

Note regarding SArKS input: the scores input to SArKS are not required to be positive, nor
are they required to be integers. Similarly, the sequences input to SArKS do not have to use
the DNA alphabet.


### 1.2.1 SArKS run parameters

Aside from input of data in the form of sequences with matching scores, sarks also requires
specification of run parameters halfWindow, spatialLength, and minGini.
I would recommend that users try several (2-4) values of halfWindow and, if interested in spa-
tial smoothing, spatialLength, using the method described in 1.6 below. For halfWindow,
you might start with halfWindow values on the order of

```
                  halfWindow ∈ {n/20, n/10, n/5}                          (1)
```

where n is the number of input sequences.

To get a sense of the intuition behind the values in (1): If we take halfWindow=n/20, we get a
full smoothing window of size 2*halfWindow+1 ≈ 10 , which amounts to looking for motifs
which would be expected to occur about once in every 10 input sequences.

The user should of course feel free to vary these fractions in either direction based on his
or her own data and goals! Especially if the number n of input sequences is very high, so
that n/20 is very large (say, in the thousands or more), these halfWindow values may be larger
than optimal.

Another consideration to keep in mind when choosing halfWindow values is that there is
a link between n, the average sequence length |w|, halfWindow, the size (really entropy)
of the sequence alphabet and the length k of the k-mer motif sequences sarks is likely
to detect. Assuming an approximately uniformly distributed alphabet (often not a valid
assumption in practice) of a distinct characters:

```
                  k ≈ log_a (n * |w| / halfWindow)                        (2)
```

should give a very rough sense of what length the motifs returned are likely to have (though
when employing spatial smoothing with merging of consecutive motifs, this will be less
accurate). Equation (2) is mosly useful for understanding how changes to halfWindow are
likely to impact output k-mer lengths, not truly predicting the exact lengths which will be
seen.

Aside from the trivial value of spatialLength=0 used to turn off spatial smoothing entirely,
it is a bit more difficult to provide general guidelines for spatialLength. Small values (say,
3–10) can be useful to help detect individual motifs even when no larger spatial structure
is really expected; in Wylie et al. [2019] we also tried spatialLength=100 when studying
regulation of gene expression as that is on the order of the low end of the length DNA
enhancer regions. Users should feel free to experiment with spatialLength values that
make sense in their own applications.

Finally, for minGini, my current advice would be to start with the value minGini=1.1 (or
minGini=0 to see what happens without this filter). Interested users should consult sections
1.4.1 and 2.2 (especially 2.2.1) as well as Wylie et al. [2019] to get a better understanding of
what this parameter does in order to get a sense of how they might choose better values for
their own applications.


## 1.3 Output

sarks then aims to identify short sequence motifs, and potentially longer sequence domains
enriched in such motifs (MMDs), where the occurrence of the identified motifs in the input
sequences is correlated with the numeric scores assigned to the sequences.
Let’s consider a minimal sarks workflow (which we’ll discuss in more detail below) to
illustrate the way in which this output is constructed):

```R
sarks <- Sarks(simulatedSeqs, simulatedScores, halfWindow=4)
filters <- sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
permDist <- permutationDistribution(sarks, reps=250, filters, seed=123)
thresholds <- permutationThresholds(filters, permDist, nSigma=2.0)
peaks <- kmerPeaks(sarks, filters, thresholds)
peaks[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*

|    i |    s | block |  wi | kmer        | windowed |
|------|------|-------|-----|-------------|----------|
| 1458 | 4442 |    24 | 175 |  ATACTGAG   |        1 |
| 1459 | 3545 |    21 |  31 |  ATACTGAGA  |        1 |
| 1460 | 3960 |    22 | 195 |  ATACTGAGA  |        1 |
| 1461 | 4519 |    25 |   1 |  ATACTGAGA  |        1 |
| 1462 | 3456 |    20 | 193 |  ATACTGAGA  |        1 |
| 1463 | 5595 |    29 |  73 |  ATACTGAG   |        1 |
| 2256 | 3544 |    21 |  30 | CATACTGAGA  |        1 |
| 2257 | 3959 |    22 | 194 | CATACTGAGA  |        1 |
| 2258 | 4518 |    25 |   0 | CATACTGAGA  |        1 |
| 5862 | 4443 |    24 | 176 |   TACTGAG   |        1 |
| 5863 | 3546 |    21 |  32 |   TACTGAGA  |        1 |
| 5864 | 3961 |    22 | 196 |   TACTGAGA  |        1 |
| 5865 | 4520 |    25 |   2 |   TACTGAGA  |        1 |


### 1.3.1 k-mers

For now let’s focus specifically on the column peaks$kmer, which tells us that SArKS has
identified

```R
unique(peaks$kmer)
```
*out:*
```
[1] "ATACTGAG"      "ATACTGAGA"     "CATACTGAGA" "TACTGAG"         "TACTGAGA"
```

as k-mers whose occurrence in simulatedSeq is associated with higher values of simulat-
edScores. sarks provides a convenience function kmerCounts:

```R
kmerCounts(unique(peaks$kmer), simulatedSeqs)
```
*out:*

|    | ATACTGAG | ATACTGAGA | CATACTGAGA | TACTGAG | TACTGAGA |
|----|----------|-----------|------------|---------|----------|
|  0 |        0 |         0 |          0 |       0 |        0 |
|  1 |        0 |         0 |          0 |       0 |        0 |
|  2 |        0 |         0 |          0 |       0 |        0 |
|  3 |        0 |         0 |          0 |       0 |        0 |
|  4 |        0 |         0 |          0 |       0 |        0 |
|  5 |        0 |         0 |          0 |       0 |        0 |
|  6 |        0 |         0 |          0 |       0 |        0 |
|  7 |        0 |         0 |          0 |       0 |        0 |
|  8 |        0 |         0 |          0 |       0 |        0 |
|  9 |        0 |         0 |          0 |       0 |        0 |
| 10 |        0 |         0 |          0 |       0 |        0 |
| 11 |        0 |         0 |          0 |       0 |        0 |
| 12 |        0 |         0 |          0 |       0 |        0 |
| 13 |        0 |         0 |          0 |       0 |        0 |
| 14 |        0 |         0 |          0 |       0 |        0 |
| 15 |        0 |         0 |          0 |       0 |        0 |
| 16 |        0 |         0 |          0 |       0 |        0 |
| 17 |        0 |         0 |          0 |       0 |        0 |
| 18 |        0 |         0 |          0 |       0 |        0 |
| 19 |        0 |         0 |          0 |       0 |        0 |
| 20 |        1 |         1 |          1 |       1 |        1 |
| 21 |        1 |         1 |          1 |       1 |        1 |
| 22 |        1 |         1 |          1 |       1 |        1 |
| 23 |        1 |         1 |          1 |       1 |        1 |
| 24 |        1 |         1 |          1 |       1 |        1 |
| 25 |        1 |         1 |          1 |       1 |        1 |
| 26 |        1 |         1 |          1 |       1 |        1 |
| 27 |        1 |         1 |          1 |       1 |        1 |
| 28 |        1 |         1 |          1 |       1 |        1 |
| 29 |        1 |         1 |          1 |       1 |        1 |

which shows us that each of the identified k-mers appears exactly once in each of the higher-
scoring sequences ("20"-"29") and never in any of the lower-scoring sequences. Looking a
bit more closely at the specific k-mers identified here, you can see that they are all
substrings of the longest one, the 10-mer "CATACTGAGA", so that the perfect agreement between
the columns of the kmerCounts output above are perhaps to be expected. sarks provides
functionality to help to identify such situations and simplify the output, as discussed further
in section 1.3.2 below.

First let’s go back to the output peaks from the call to kmerPeaks (focusing on the 8th row
of peaks):

```R
peak8 = peaks[8, ]
peak8[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*

|    i |    s | block |  wi | kmer       | windowed |
|------|------|-------|-----|------------|----------|
| 2257 | 3959 |    22 | 194 | CATACTGAGA |        1 |

and consider the columns:

- **i**  
  suffix array index
- **s**  
  suffix array value si
- **block**  
  the name of the input sequence from which the indicated k-mer is derived
- **wi**  
  the (0-based!) position ωi of the k-mer within the indicated block/sequence

Considering block and wi first, this tells us that, setting

```R
block <- simulatedSeqs[[peak8$block]]
kmerStart <- peak8$wi + 1
kmerEnd <- kmerStart + nchar(peak8$kmer) - 1
```

(noting the + 1 required to define kmerStart in 1-based R relative to 0-based wi) we should
find that

```R
substr(block, kmerStart, kmerEnd)
```
*out:*
```
[1] "CATACTGAGA"
```

yields the same value as
```R
peak8$kmer
```
*out:*
```
[1] "CATACTGAGA"
```


### 1.3.2 Reducing redundancy in peaks

In cases such as "CATACTGAGA" in the simulated data set analyzed here, sarks can simplify
k-mer output using a couple of auxiliary functions:

```R
nonRedundantPeaks <- mergedKmerSubPeaks(sarks, filters, thresholds)
nonRedundantPeaks[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*

|    i |    s | block |  wi | kmer       | windowed |
|------|------|-------|-----|------------|----------|
| 1462 | 3456 |    20 | 193 |  ATACTGAGA |        1 |
| 2256 | 3544 |    21 |  30 | CATACTGAGA |        1 |
| 2257 | 3959 |    22 | 194 | CATACTGAGA |        1 |
| 1458 | 4442 |    24 | 175 |  ATACTGAG  |        1 |
| 2258 | 4518 |    25 |   0 | CATACTGAGA |        1 |
| 1463 | 5595 |    29 |  73 |  ATACTGAG  |        1 |

Here mergedKmerSubPeaks removes redundant k-mer output where multiple k-mer peaks
are reported with successive spatial s coordinates (e.g., the reported peak peaks[3, ]
from block="22" at wi=195 and kmer="ATACTGAGA" immediately follows peaks[8, ] with
block="22" and wi=194 and kmer="CATACTGAGA", and is thus merged with that peak).
Further simplifiation is still possible in this case using:

```R
extendedPeaks <- extendKmers(sarks, nonRedundantPeaks)
extendedPeaks[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*

|    i |    s | block |  wi | kmer       | windowed |
|------|------|-------|-----|------------|----------|
| 2259 | 3455 |    20 | 192 | CATACTGAGA | NA       |
| 2256 | 3544 |    21 |  30 | CATACTGAGA | 1        |
| 2257 | 3959 |    22 | 194 | CATACTGAGA | 1        |
| 2255 | 4441 |    24 | 174 | CATACTGAGA | NA       |
| 2258 | 4518 |    25 |   0 | CATACTGAGA | 1        |
| 2260 | 5594 |    29 |  72 | CATACTGAGA | NA       |

which detects that even though the full occurrence of "CATACTGAGA" in input sequence "24"
was not reported as a k-mer by SArKS, the k-mer "ACACTGAG" in nonRedundantPeaks[1, ] is
flanked in sequence "24" by a "C" to the left and an "A" to the right and can hence be
extended to math the full motif.


### 1.3.3 Suffix arrays

How are we to interpret the columns i and s in peaks? This requires understanding a bit
more about SArKS: Specifically, that the method begins by concatenating all of the input
sequences into one big string (called catSeq in the underlying java object referenced by
sarks):

```R
concatenated <- sarks$getCatSeq()
nchar(concatenated)
```
*out:*
```
[1] 7530
```

concatenated is actually a bit longer than the sum of the lengths of the input sequences
because it keeps track of where one sequence ends and another begins using a special
(dollar-sign) character. In this way, concatenated is divided into separate “blocks,” each
corresponding to one of the input sequences.

The column s in peaks then gives us the (0-based!) position of the annotated k-mer in the
concatenated sequence:

```R
kmerCatStart <- peak8$s + 1
kmerCatEnd <- kmerCatStart + nchar(peak8$kmer) - 1
substr(concatenated, kmerCatStart, kmerCatEnd)
```
*out:*
```
[1] "CATACTGAGA"
```

But what about i? As we will confirm, it is the (0-based) position of the suffix

```R
theSuffix <- substr(concatenated, kmerCatStart, nchar(concatenated))
```
*out:*
in the list of all suffixes after they have been lexicographically sorted :

```R
extractSuffix <- function(s) {
    ## returns suffix of concatenated starting at position s (1-based)
    substr(concatenated, s, nchar(concatenated))
}
allSuffixes <- vapply(1:nchar(concatenated), extractSuffix, '')
sortedSuffixes <- sort(allSuffixes)
```

Sure enough,

```R
i1based <- peak8$i + 1
sortedSuffixes[i1based] == theSuffix
```
*out:*
```
[1] TRUE
```


### 1.3.4 Window averaging (kernel smoothing)

Because the suffixes in sortedSuffixes are sorted, the suffixes in a window centered on
i1based all start with the same few characters (a k-mer):

```R
iCenteredWindow <- (i1based - 4):(i1based + 4)
iCenteredWindowSuffixes <- sortedSuffixes[iCenteredWindow]
all(substr(iCenteredWindowSuffixes, 1, 10) == 'CATACTGAGA')
```
*out:*
```
[1] TRUE
```

For each suffix in sortedSuffixes, we can identify which input sequence contributed the
block of concatenated where the suffix begins:

```R
iCenteredWindow0Based <- iCenteredWindow - 1
sourceBlock(sarks, i=iCenteredWindow0Based)
```
*out:*
```
[1] "26" "28" "24" "21" "22" "25" "20" "29" "27"
```

Note that the 9 suffixes in iCenteredWindowSuffixes derive from come from 9 of the 10
higher-scoring sequences (all but "23"). This is not an accident: since the motif "CATACT-
GAGA" is only present in high-scoring sequences, the suffixes starting with "CATACTGAGA"
must derive from blocks associated with high-scoring sequences. Thus the average score
of the sequences contributing the suffixes specified by iCenteredWindowSuffixes must be
high (value of 1) as well.

The SArKS algorithm turns this around to identify motifs by hunting for windows in the
sorted suffix list where the average score of the corresponding sourceBlock sequences is
high.

The windowed column of the peaks output contains the average sourceBlock sequence
scores for the windows centered around each sorted suffix position i and extending to both
the left and right by halfWindow=4 positions. These average values are also referred to as
ŷi in Wylie et al. [2019]; a vector containing all of these values can be obtained from the
object sarks:

```R
yhat <- sarks$getYhat()
i0based <- seq(0, length(yhat)-1)
plot(i0based, yhat, type='l', xlab='i')
```
*output plot omitted; see pdf version of vignette in R package*

From the plot of ŷi against i, we can see three peaks at which ŷi spikes up to 1,
corresponding to the ranges i ∈ [1458, 1463], i ∈ [2256, 2258], and i ∈ [5862, 5865]
which make up the values in peaks$i.


## 1.4 Permutations and thresholds

How do we know that a given value of ŷi (or, equivalently, peaks$windowed) is high enough to
include in the peaks output? In order to set such a threshold without having to make possibly
unwarranted assumptions about the structure of the input sequences or the distribution of
the scores, SArKS employs a permutational approach.
To illustrate the idea here, let’s randomly permute the scores—so that the permuted scores
have no relationship with which sequences contain "CATACTGAGA"—and construct a new
Sarks object using the permuted scores:

```R
set.seed(12345)
scoresPerm <- sample(simulatedScores)
names(scoresPerm) <- names(simulatedScores)
sarksPerm <- Sarks(simulatedSeqs, scoresPerm, halfWindow=4)
```

If we try using the same procedure we did before to look for k-mer peaks:

```R
permDistNull <- permutationDistribution(
        sarksPerm, reps=250, filters, seed=123)
thresholdsNull <- permutationThresholds(
        filters, permDistNull, nSigma=2.0)
peaksNull <- kmerPeaks(sarksPerm, filters, thresholdsNull)
peaksNull[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*
```
[1] i        s        block    wi                   kmer          windowed
<0 rows> (or 0-length row.names)
```

we find that SArKS does not detect anything worthy of reporting after we disrupt any
association between input sequences and input scores by permuting the scores.

This is to be expected given that the thresholds θ used by SArKS are defined by considering
many (here, reps=250) such permutations and choosing θ such that only very few random
permutations would produce smoothed ŷi scores exceeding θ. The adjustable parameter
nSigma controls how stringent the thresholds are: higher nSigma values lead to higher
thresholds, reducing sensitivity but also reducing the false positive rate.

Once we’ve set thresholds using permutationThresholds, we can estimate the false positive
rate, defined here as the frequency of seeing nonempty k-mer result sets when there is really
no association between sequence and score:

```R
fpr = estimateFalsePositiveRate(
        sarks, reps=250, filters, thresholds, seed=321)
fpr$ci
```
*out:*
```
    method x   n mean      lower      upper
1    exact 5 250 0.02 0.00652507 0.04605358
```

indicating a point estimate of 2% for the false positive rate, with a 95% confidence interval
of (0.65%, 4.6%). This confidence interval can be made tighter by using a higher number for
reps in the estimateFalsePositiveRate call.

Note regarding random number generator seed for estimateFalsePositiveRate: do not use
the same seed for estimateFalsePositiveRate as was used in permutationDistribution
call used to set thresholds. The false positive rate estimation procedure assumes that the
random permutations used in estimateFalsePositiveRate are independent of those used
to set thresholds.


### 1.4.1 Gini impurity filter

When considering how to set SArKS thresholds in order to obtain a reasonable tradeoff
between sensitivity and false positive rate, the mysterious minGini parameter should be
factored in as well. This parameter is discussed in more detail in sections 2.1.7 and 2.2.1, but
the basic idea is to filter likely false positive suffix array index positions i out of
consideration regardless of their smoothed scores ŷi . These likely false positive indices i
come from excessive repetition of the same input sequences contributing repeatedly to the
smoothing windows centered on i.

My recommendation is to set minGini to a value slightly greater than 1: the value 1.1 seems
empirically to work well in many situations. Note that for minGini > 1, this filter becomes
more stringent the closer to 1 it is set; the opposite is true if you set minGini to a value
less than 1, and you can turn the filter off entirely by setting it to 0. See section 2.2.1
and Wylie et al. [2019] for more details.

Here we examine the effects of this parameter using the simulated data; first with the filter:

```R
estimateFalsePositiveRate(
    sarks, reps=250, filters, thresholds, seed=321)$ci
```
*out:*
```
    method x   n mean      lower      upper
1    exact 5 250 0.02 0.00652507 0.04605358
```

And now without:

```R
filtersNoGini <- sarksFilters(halfWindow=4, spatialLength=0, minGini=0)
estimateFalsePositiveRate(
    sarks, reps=250, filtersNoGini, thresholds, seed=321)$ci
```
*out:*
```
    method x    n mean      lower     upper
1    exact 42 250 0.168 0.1238429 0.2202254
```

Of course we could compute a new set of thresholdsNoGini using filtersNoGini, but the
results of this on our ability to detect anything are worth considering:

```R
permDistNoGini <-
        permutationDistribution(sarks, reps=250, filtersNoGini, seed=123)
thresholdsNoGini <-
        permutationThresholds(filtersNoGini, permDistNoGini, nSigma=2.0)
peaksNoGini <- kmerPeaks(sarks, filtersNoGini, thresholdsNoGini)
peaksNoGini[ , c('i', 's', 'block', 'wi', 'kmer', 'windowed')]
```
*out:*
```
[1] i        s        block    wi                kmer     windowed
<0 rows> (or 0-length row.names)
```

### 1.4.2 Multithreading sarks permutational analyses

Setting thresholds and estimating false positive rates using permutation methods is usually
the most time-consuming step in the SArKS workflow. To facilitate more rapid turnaround,
the java back end of sarks supports multithreading for these processes. In order to take
advantage of multithreading, you just need to specify how many threads you’d like to use
when you invoke the Sarks constructor using the nThreads argument to that function; all
permutation steps performed using the resulting sarks object will then use the specified
number of threads.

Note regarding sarks multithreading: while the different threads performing the permuta-
tional analyses share as many data structures as possible to reduce the memory requirements
of sarks, each thread will need to keep track of its own permuted smoothed scores (and
spatially smoothed scores if necessary). This can increase memory requirements quickly, so
make sure you are running sarks on a machine with large RAM if you are planning to use
many threads on a large data set.


## 1.5 Spatial smoothing

So far we have not taken advantage of the spatial smoothing features in SArKS. While one
of the major uses of spatial smoothing is to detect multi-motif domains (MMDs), which
are not present in the simple simulated data set included in the sarks package, spatial
smoothing can also help sharpen the ability of SArKS to detect individual motifs:

```R
sarks <- Sarks(
    simulatedSeqs, simulatedScores, halfWindow=4, spatialLength=3
)
filters <- sarksFilters(halfWindow=4, spatialLength=3, minGini=1.1)
permDist <- permutationDistribution(sarks, reps=250, filters, seed=123)
thresholds <- permutationThresholds(filters, permDist, nSigma=5.0)
## note setting of nSigma to higher value of 5.0 here!
peaks <- kmerPeaks(sarks, filters, thresholds)
peaks[ , c('i', 's', 'block', 'wi', 'kmer', 'spatialWindowed')]
```
*out:*

|    i |    s | block |  wi | kmer       | spatialWindowed |
|------|------|-------|-----|------------|-----------------|
| 1458 | 4442 |    24 | 175 |  ATACTGAG  |       0.9259259 |
| 1459 | 3545 |    21 |  31 |  ATACTGAGA |       0.9629630 |
| 1460 | 3960 |    22 | 195 |  ATACTGAGA |       0.9259259 |
| 1461 | 4519 |    25 |   1 |  ATACTGAGA |       0.9259259 |
| 2256 | 3544 |    21 |  30 | CATACTGAGA |       1.0000000 |
| 2257 | 3959 |    22 | 194 | CATACTGAGA |       1.0000000 |
| 2258 | 4518 |    25 |   0 | CATACTGAGA |       1.0000000 |
| 5863 | 3546 |    21 |  32 |   TACTGAGA |       0.9259259 |

Note that here we use nSigma=5.0, a much more stringent setting than the nSigma=2.0
value used when no spatial smoothing was employed; had we tried nSigma=5.0 without
spatial smoothing, SArKS would not have been able to detect the "CATACTGAGA" motif.

The peaks object returned by kmerPeaks when spatial smoothing is employed contains the i
and s coordinates for the left endpoints of spatial windows (windows in s-space, not i-space)
enriched in high ŷ scores. Especially when these windows have longer spatialLength values,
it can also be useful to pick out individual motif “subpeaks” within these spatial windows:

```R
subpeaks <- mergedKmerSubPeaks(sarks, filters, thresholds)
subpeaks[ , c('i', 's', 'block', 'wi', 'kmer')]
```
*out:*

|    i |    s | block |  wi | kmer       |
|------|------|-------|-----|------------|
| 2256 | 3544 |    21 |  30 | CATACTGAGA |
| 2257 | 3959 |    22 | 194 | CATACTGAGA |
| 1458 | 4442 |    24 | 175 |  ATACTGAG  |
| 2258 | 4518 |    25 |   0 | CATACTGAGA |

In this case, with spatialLength=3, this is not such an important step—although it does
serve to accomplish some simplification of results as was described when mergedKmerSub-
Peaks was first introduced in section 1.3.2—but in general it is a critical piece of the SArKS
methodology when spatial smoothing is in use.

The subpeaks output of mergedKmerSubPeaks should generally be regarded as the main
individual motif output for SArKS. Section 2.5 provides more details.
When employing spatial smoothing to identify MMDs, you may want to inspect individual
input sequences of interest by visualizing the SArKS results:

```R
block22 <- blockInfo(sarks, block='22', filters, thresholds)
library(ggplot2)
ggo <- ggplot(block22, aes(x=wi+1)) ## +1 because R indexing is 1-based
ggo <- ggo + geom_point(aes(y=windowed), alpha=0.6, size=1)
ggo <- ggo + geom_line(aes(y=spatialWindowed), alpha=0.6)
ggo <- ggo + geom_hline(aes(yintercept=spatialTheta), color='red')
ggo <- ggo + ylab('yhat') + ggtitle('Input Sequence "22"')
ggo <- ggo + theme_bw()
print(ggo)
```
*output plot omitted; see pdf version of vignette in R package*

which here clearly shows the spike at the location wi+1=195 of the "CATACTGAGA" motif in
sequence "22".

The use of the more stringent nSigma=5.0 setting reduces the false positive rate relative to
the nSigma=2.0 setting we used when not employing spatial smoothing:

```R
estimateFalsePositiveRate(
    sarks, reps=250, filters, thresholds, seed=321)$ci
```
*out:*
```
    method x   n mean lower      upper
1    exact 0 250    0     0 0.01464719
```


## 1.6 Varying SArKS parameters

You may be wondering what the point of the filters object created in the SArKS workflow
is, as it seems to specify the halfWindow ad spatialLength parameters in a manner
redundant to their specification in the Sarks constructor. Aside from the specification of the
minGini parameter, the motivation for this step is that we usually don’t know exactly what
halfWindow or spatialLength values will yield the most useful output (and the answer to
these questions will be different for different data sets and different questions).

To address this, we can use sarksFilters to test a variety of combinations of halfWindow
and spatialLength (and, if so desired, minGini as well) values like so:

```R
filters <- sarksFilters(
        halfWindow=c(4, 8), spatialLength=c(2, 3), minGini=1.1)
permDist <- permutationDistribution(sarks, reps=250, filters, seed=123)
thresholds <- permutationThresholds(filters, permDist, nSigma=5.0)
peaks <- mergedKmerSubPeaks(sarks, filters, thresholds)
peaks[ , c('halfWindow', 'spatialLength', 'i', 's', 'block', 'wi', 'kmer')]
```
*out:*

| halfWindow | spatialLength |    i |    s | block |  wi | kmer       |
|------------|---------------|------|------|-------|-----|------------|
|          4 |             2 | 2256 | 3544 |    21 |  30 | CATACTGAGA |
|          4 |             2 | 2257 | 3959 |    22 | 194 | CATACTGAGA |
|          4 |             2 | 1458 | 4442 |    24 | 175 |  ATACTGAG  |
|          4 |             2 | 2258 | 4518 |    25 |   0 | CATACTGAGA |
|          4 |             3 | 2256 | 3544 |    21 |  30 | CATACTGAGA |
|          4 |             3 | 2257 | 3959 |    22 | 194 | CATACTGAGA |
|          4 |             3 | 1458 | 4442 |    24 | 175 |  ATACTGAG  |
|          4 |             3 | 2258 | 4518 |    25 |   0 | CATACTGAGA |
|          8 |             2 | 1459 | 3545 |    21 |  31 |  ATACTGA   |
|          8 |             2 | 2261 | 5173 |    27 | 153 | CATACTGA   |
|          8 |             3 | 5860 | 4947 |    26 | 178 |   TACTGA   |
|          8 |             3 | 2261 | 5173 |    27 | 153 | CATACTGA   |

Note that, as is generally true with SArKS, the results obtained using larger halfWindow
values tend to have shorter identified k-mers.

As always, testing multiple parameter sets can increase false positive rates:

```R
estimateFalsePositiveRate(
    sarks, reps=250, filters, thresholds, seed=321)$ci
```
*out:*
```
    method x   n mean         lower      upper
1    exact 1 250 0.004 0.0001012661 0.02208387
```

As suggested in section 1.5 above, this can be countered by using increased nSigma values
when setting thresholds if we plan to test a wider range of possible SArKS parameters.
Note regarding multiple hypothesis testing in SArKS: The false positive rates estimated by
SArKS test the rate of detecting any results using any of the specified parameter settings, and
are thus a type of family-wise error rate. As long as all parameter combinations tested are
included in the filters employed in estimateFalsePositiveRate, no further adjustment
for multiple testing should be applied to the estimated false positive rates.


## 1.7 Clustering similar k-mers into broader motifs

While there is essentially one unambiguous k-mer "CATACTGAGA" of note in the simulated-
Seqs-simulatedScores example, real data sets will generally have more complex sequence
motifs allowing for some variation both in length and composition of the patterns. When
applying SArKS methodology to real data, many similar k-mers representing the same basic
motif may result: sarks provides functions for clustering these k-mers into broader motif
patterns.

For example, let’s say we had obtained the following k-mers using sarks:

```R
kmers <- c(
    'CAGCCTGG', 'CCTGGAA', 'CAGCCTG', 'CCTGGAAC', 'CTGGAACT',
    'ACCTGC', 'CACCTGC', 'TGGCCTG', 'CACCTG', 'TCCAGC',
    'CTGGAAC', 'CACCTGG', 'CTGGTCTA', 'GTCCTG', 'CTGGAAG', 'TTCCAGC'
)
```

We could cluster these like so:

```R
kmClust <- clusterKmers(kmers, directional=FALSE)
## directional=FALSE indicates that we want to treat each kmer
##                   as equivalent to its reverse-complement
kmClust
```
*out:*
```
$CAGCCTGG
[1] "CAGCCTGG" "CAGCCTG"

$CTGGAAC
[1] "CCTGGAA"     "CCTGGAAC" "CTGGAACT" "TCCAGC"           "CTGGAAC"    "CTGGAAG"     "TTCCAGC"

$CACCTGC
[1] "ACCTGC"     "CACCTGC" "CACCTG"      "CACCTGG"

$TGGCCTG
[1] "TGGCCTG"

$CTGGTCTA
[1] "CTGGTCTA"

$GTCCTG
[1] "GTCCTG"
```

The resulting object kmClust is a named list: each element of this list is a character vector
listing the elements of the vector kmers composing the corresponding cluster, while the name
of the cluster is a k-mer from kmers found to be particularly representative of the cluster.
sarks then allows us to count how many times each motif (=cluster of k-mers) occurs in
each sequence:

```R
clCounts <- clusterCounts(kmClust, simulatedSeqs, directional=FALSE)
## directional=FALSE to count hits of either a kmer from the cluster
##                   or the reverse-complement of such a kmer
## clCounts is a matrix with:
## - one row for each sequence in simulatedSeqs
## - one column for each *cluster* in kmClust
head(clCounts)
```
*out:*

| CAGCCTGG | CTGGAAC | CACCTGC | TGGCCTG | CTGGTCTA | GTCCTG |
|----------|---------|---------|---------|----------|--------|
|        0 |       0 |       0 |       0 |        0 |      0 |
|        0 |       0 |       1 |       0 |        0 |      0 |
|        0 |       0 |       0 |       0 |        0 |      0 |
|        0 |       0 |       0 |       0 |        0 |      0 |
|        0 |       0 |       0 |       0 |        0 |      0 |
|        0 |       0 |       0 |       1 |        0 |      1 |

Can also get specific information about the location of these matches:

```R
clLoci <- locateClusters(
    kmClust, simulatedSeqs, directional=FALSE, showMatch=TRUE
)
## showMatch=TRUE includes column specifying exactly what k-mer
##                registered as a hit;
##                this can be very slow, so default is showMatch=FALSE
clLoci
```
*out:*

| seqid | cluster  | location | match    |
|-------|----------|----------|----------|
|    25 | CAGCCTGG |      143 | CAGCCTGG |
|     6 | CTGGAAC  |      190 | GCTGGAA  |
|    16 | CTGGAAC  |      196 | AGTTCCAG |
|    28 | CTGGAAC  |       36 | GTTCCAG  |
|     1 | CACCTGC  |      242 | GCAGGTG  |
|    28 | CACCTGC  |       25 | GCAGGT   |
|     5 | TGGCCTG  |      216 | TGGCCTG  |
|    19 | CTGGTCTA |       51 | TAGACCAG |
|     5 | GTCCTG   |       19 | GTCCTG   |
|     9 | GTCCTG   |       67 | GTCCTG   |
|    13 | GTCCTG   |       86 | CAGGAC   |
|    15 | GTCCTG   |       19 | CAGGAC   |
|    17 | GTCCTG   |      151 | CAGGAC   |
|    19 | GTCCTG   |       76 | CAGGAC   |
|    24 | GTCCTG   |       38 | GTCCTG   |

This shows us that, for example, there is one match for the cluster "CAGCCTGG" spanning
sequence characters 143-150 of sequence "25", and that this is an exact match of the k-mer
"CAGCCTGG" in its forward orientation.

These results also tells us that there are three matches of k-mers from the cluster "CTGGAAC",
in sequences "6", "16", and "28". The specific k-mers found are different in each case:
- sequence "6" has a hit for the reverse-complement of k-mer "TTCCAGC", while
- sequence "16" has a hit for k-mer "CTGGAACT" in reverse-complement orientation and
- sequence "28" has a hit for "CTGGAAC" also in reverse-complement orientation.


# 2 How sarks works

*Omitted because includes extensive equations;
see pdf vignette included in R package for this section*


# 3 Notation glossary

*Omitted because includes extensive equations;
see pdf vignette included in R package for this section*


# References
- J. Kärkkäinen and P. Sanders. Simple linear work suffix array construction. In International
  Colloquium on Automata, Languages, and Programming, pages 943–955. Springer, 2003.

- D. Wylie, H. Hofmann, and B. Zemelman. SArKS: de novo discovery of gene expression regulatory
  motif sites and domains by suffix array kernel smoothing. Bioinformatics, 35(20):3944–3952, 2019.
