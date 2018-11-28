# sarks user-guide

Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
multi-motif domains (MMDs)
(https://www.biorxiv.org/content/early/2018/10/25/133934) whose
presence or absence within each of a set of sequences (contained
within a fasta file) is correlated with numeric scores assigned to
those sequences (contained with a two column tab-separated-values
file, or tsv).

For installation instructions, consult [INSTALL.md](INSTALL.md).

# sarks command-line execution

Included in the sarks project examples/ folder are several scripts for
usage at the command line:

- [**sarkselect.py**](#sarkselectpy)  
  main script for motif discovery
- [**sarkstest.py**](#sarkstestpy)  
  estimate false positive rate associated with motif set determined by
  sarkselect.py
- [**extract_kmers.py**](#extract\_kmerspy)  
  extract list of unique k-mers identified by sarkselect.py
- [**zrank_singly_smoothed_kmers.py**](#zrank\_singly\_smoothed\_kmerspy)  
  rank those k-mers which were selected without use of spatial
  smoothing according to maximum z threshold at which they could be
  identified
- [**cluster_seqs.R**](#cluster\_seqsr)  
  assign k-mers to clusters; requires specification of number of
  clusters to generate
- [**count_kmers.py**](#count\_kmerspy)  
  count (or locate) occurrences of k-mer motifs or of clusters of
  k-mers within sequences in a fasta file

All Python scripts are written for Python 3. Details regarding the
usage of these scripts are provided below (abbreviated usage
information can be obtained for each script by calling it with the -h
or --help flag).

Examples of the usage of these scripts can be found in the files
indicated below (I suggest going through them in the indicated order,
reading through the comments and running the commands line-by-line):

- [examples/simulated_example.md](examples/simulated_example.md)  
  toy simulated data example to illustrate basic functionality of
  sarks
- [examples/mo2015_downstream_example.md](examples/mo2015_downstream_example.md)  
  minimal introductory example of sarks applied to RNA-seq
  differential expression results, can still be run fairly quickly
  (though final sarkstest.py step may take 10+ minutes depending on
  hardware)
- [examples/mo2015_upstream_example.md](examples/mo2015_upstream_example.md)  
  this example uses longer sequences and includes spatial smoothing,
  hence it takes longer (hours) and requires more RAM (**recommend
  minimum 48G available**)

## sarkselect.py

Main input for sarkselect.py is provided by the first two arguments
listed below: a fasta file containing the sequences and a score file
containing the scores for those sequences.

- **-f**  
  fasta file containing sequences to analyze
- **-s**  
  scores tab-separated values (tsv) file:
    - column 1: sequence identifiers matching the fasta file
    - column 2: numeric scores
- **-w**  
  half window width for first kernel smoothing pass (kappa in paper),
  can supply multiple values using commas (no spaces)
- **-l**  
  spatial smoothing length (lambda in paper),
  can supply multiple values using commas (no spaces)  
  default: 0 (no spatial smoothing)
- **-g**  
  parameter for calculation of Gini impurity filter (gamma in paper),
  can supply multiple values using commas (no spaces)  
  default: 1.1
- **-r**  
  number R of permutations to use in setting significance thresholds
  for peak-calling  
  default: 100
- **-z**  
  multiple z of standard deviations above mean (of maximum smoothed
  suffix scores obtained after randomly permuting scores assigned to
  sequences) defining threshold (see supplementary Section S2.6, Eq
  (S24-S25) of paper)  
  default: 4.0
- **-o**  
  output directory to be created/overwritten

Additional optional flags/arguments:

- **-p**  
  takes no argument;
  logical flag indicating to remove redundant kmer output,
  as described in supplementary Section S2.3 of the paper.
  Most useful for small data sets  
  default: False
- **-n**  
  takes no argument; logical flag indicating to find motifs associated
  with lowest/most negative scores instead of highest/most positive scores  
  default: False
- **-c**  
  specifies path for outputting concatenated fasta file  
  default: a temporary file is automatically created (and removed when done)
- **-e**
  fixes random number generator seed to generate reproducible output  
  default: no seed

## sarkstest.py

sarkstest.py requires output from sarkselect.py.

Most of the arguments for sarkstest.py are identical to those of
sarkselect.py: sarkstest.py should be called (at least initially, see
more below) using the same set of arguments as sarkselect.py was
*except* for the -r (reps) and the -o (output) flag:

- **-r**  
  number R_2 of permutations used for false positive rate testing,
  need not be same as number R of permutations used to set threshold  
  default: 100
- **-i**  
  input directory for sarkstest.py is output directory from sarkselect.py;
  needed to determine selection thresholds for testing
- **-o**  
  output directory for sarkstest.py; will be created/overwritten

Aside from the output saved to the specified output directory,
sarkstest.py also emits three lines of text on stdout:
1. how many permutations led to non-empty motif sets using the
   thresholds determined by the sarkselect.py run;
2. resulting point estimate of false positive rate;
3. estimated 95% confidence interval (CI) for false positive rate.

After sarkstest.py has been run once on the output from sarkselect.py,
you can re-run just to get these false positive rate estimates without
recomputing permutations once they've been generated. In this case,
only the -z, -i, and -o flags should be passed, with -o set to the
same directory as was previously output by sarkstest.py.

## extract_kmers.py

Takes one positional argument: tab-separated values (tsv) file with a
header row in which one column is named "kmer" (no quotation marks
required in the column name). The files peaks.tsv or merged_peaks.tsv
(if spatial smoothing was employed) produced by sarkselect.py (in the
resulting output directory) are good candidates.

Outputs (to stdout) a list of unique k-mers (one per line) found in
the kmer column of the input file. "Unique" is here taken by
default to consider two k-mers to be equivalent if one is the
reverse complement of the other; if instead you want k-mers treated
as distince unless they are exactly equal as strings, you can pass
the **-d** (directional) flag at the command line.

**NOTE:** if the -d flag is not activated, the lexicographically lower
(alphabetically first) orientation of each k-mer will be reported,
*not necessarily the same orientation initially identified by SArKS*.

## zrank_singly_smoothed_kmers.py

Takes one positional argument, which should be path to output
directory from sarkselect.py.

One approach to ranking the set of k-mers identified by SArKS is
according to what the maximum value of the -z threshold parameter in
sarkselect.py at which they would be identified.

This approach is much simpler to implement for those k-mers which were
selected without spatial smoothing: When spatial smoothing is
employed, the process for merging adjacent peaks makes disrupts this
method of motif ranking.

This script thus considers only those k-mers which were selected without
spatial smoothing (hence "singly_smoothed"), returning (on stdout)
a tab-delimited table with columns

- **kmer**  
- **halfWindow**  
  value of the smoothing half-window (-w) for which k-mer was selected
- **minGini**  
  value of the gamma (-g) parameter used to select Gini impurity
  threshold for which k-mer was selected
- **zmax**  
  the maximum value of -z threshold parameter for which k-mer would
  have been identified

**NOTE:** zmax values from repeated runs of sarkselect.py will vary
somewhat because permutations used to assess null distribution are
stochastic.

## cluster_seqs.R

**NOTE:** requires R installation including the three libraries:
- **philentropy** (https://cran.r-project.org/web/packages/philentropy/index.html)
- **msa** (https://bioconductor.org/packages/release/bioc/html/msa.html)
- **cluster** (https://cran.r-project.org/web/packages/cluster/index.html)

Takes two positional arguments:

1. tab-separated values (tsv) file with a header row in which one
   column is named "kmer" (no quotation marks required in the column
   name), such as the files peaks.tsv or merged_peaks.tsv (if spatial
   smoothing was employed) produced by sarkselect.py (in the resulting
   output directory).
2. number of clusters desired.

Can add optional flag **--pfm** (no argument) to reformat cluster
output as position frequency matrices (PFMs).

Outputs a multiple sequence alignment (or PFM summary thereof) for
each cluster, with clusters separated by blank lines, to stdout.

Hierarchical clustering (average linkage) is performed based on
Jaccard coefficient distance metric applied treating each k-mer as the
set of all tetramers which can be found as substrings within it.
Tetramers which are reverse-complements of each other are treated as
equivalent, and input k-mers may be reported in reverse-complement
orientation in clustered output to facilitate multiple sequence
alignment.

For each cluster, the k-mer for which the average Jaccard distance to
all other k-mers in the cluster is minimal is selected as the
representative k-mer of the cluster and is used to name the cluster.

## count_kmers.py

Takes two files as input, one specifying either the set of k-mers or
the set of clusters of k-mers to count (or locate) and the other the
set of sequences in which to count (or locate) them.

- **-k**  
  Input file containing k-mers to count/locate. Either a tab-separated
  values (tsv) with a column named "kmer" (no quotation marks in
  column name) or a simple list of k-mers, one per line  
  (alternately, can provide comma-delimited list (no spaces) of k-mers
  directly on command line).
- **-c**  
  Alternate input file containing k-mer clusters (formatted as in the
  output from clusters_seqs.R above) to count/locate. Only one of
  -k or -c should be provided.
- **-f**  
  Input fasta file containing sequences in which k-mers are to be
  counted or located.

Returns tab-delimited tabular output on stdout: first column is sequence
identifier, then add one column for each k-mer or k-mer cluster.
  
Optional flags:
  
- **-d**  
  Specifies *directional* output: count k-mers only in the forward
  (i.e., not reverse-complemented) orientation
- **-l**  
  Species *locational* output: return tab-delimited tabular output
  with columns for
    1. **seqid**
    2. **location**
    3. **kmer** (or **cluster**)
- **-o**  
  Allow *overlapping* occurrences of k-mer input; not available for
  clustered input.

# sarks Python module

The core functionality of sarks is implemented as a Python module
defining a class *Sarks*:

```python
from sarks import Sarks
```

This class implements suffix array kernel smoothing for de novo
correlative motif discovery, conducted in several steps:
1. initialize Sarks object (requires fasta file and pandas Series
   containing scores, as well as specification of smoothing window
   size),
2. assessment of distribution of scores under null-hypothesis using
   permutationDistribution specifying Gini impurity cutoffs (minGini
   and/or minSpatialGini),
3. call to peaks method specifying filter parameters (e.g., minGini
   and smoothed score cutoff theta),
4. extraction of kmer table for top peaks using subtable method
   (followed by optional use of extendKmers if desired),
5. clustering of extracted kmers using clusterKmers.

Construction of a Sarks object requires:

- **inFasta**  
  specification of fasta file containing sequences to be analyzed
- **catFasta**  
  file name to be used for fasta file containing concatenated sequence
  produced by Sarks
- **scores**  
  pandas Series object with index matching the sequence ids in inFasta
- **halfWindow**  
  half-width of smoothing window

For example, using the example files simulated\_seqs.fa and
simulated\_scores.tsv:

```python
import pandas as pd
from sarks import Sarks

# load scores into pandas Series
scores = pd.read_table('simulated_scores.tsv', header=0, index_col=0).iloc[:, 0]
# initialize Sarks object
sarks = Sarks(
    inFasta = 'simulated_seqs.fa',          ## input sequences;
                                            ##   names should match scores.index
    catFasta = 'simulated_seqs_concat.fa',  ## name for intermediate
                                            ##   concatenated fasta file
                                            ##   (will be created/overwritten)
    scores = scores,
    halfWindow = 4                          ## kernel half-width kappa
)

```

Once one has constructed such an object, you can use its *peaks*
method to identify suffix positions corresponding to smoothed score
peaks

```python
topLocations = sarks.peaks(theta = 1,
                           minGini = 0)
```

The full list of arguments for peaks are:

- **theta**  
  minimum smoothed score value for peak to be reported (float;
  required)
- **prune**  
  whether to apply **pruneIntervals** to resulting kmer intervals
  (bool; default False)
- **spatialLength**  
  length of spatial window to apply (if any) (int; default None)
- **spatialTheta**  
  minimum spatially smoothed score value for peak to be reported (if
  any; requires spatialLength to be specified) (float; default None)
- **k**  
  kmer length to consider (int; default 12)
- **minK**  
  kmer length to report (int; default None)
- **minGini**  
  minimum Gini impurity value for suffix position to be reported as
  peak (if less than 1.0) *or* gamma value used to determine minimum
  Gini impurity value (if greater than 1.0) (float; default None)
- **minSpatialGini**  
  minimum spatially-averaged Gini impurity value for suffix position
  to be reported as peak (float; default None)
- **deduplicate**  
  if True, report only highest scoring position corresponding to any
  given kmer (bool, default False)

The locations of the peaks in and of themselves are not very
interesting; you can use the sarks objects' *subtable* method to get
more detailed information about these peaks:

```python
topTable = sarks.subtable(topLocations).sort_values('khat',
                                                    ascending=False)
```

topTable then contains

|    i |    s | kmer       |   khat | block |  wi |     gini | score | windowed |
|------|------|------------|--------|-------|-----|----------|-------|----------|
| 2257 | 3959 | CATACTGAGA | 10.250 |    22 | 194 | 0.888889 |     1 |      1.0 |
| 2258 | 4518 | CATACTGAGA | 10.250 |    25 |   0 | 0.888889 |     1 |      1.0 |
| 2256 | 3544 | CATACTGAGA |  9.625 |    21 |  30 | 0.864198 |     1 |      1.0 |
| 1460 | 3960 | ATACTGAGA  |  9.250 |    22 | 195 | 0.888889 |     1 |      1.0 |
| 1461 | 4519 | ATACTGAGA  |  9.250 |    25 |   1 | 0.888889 |     1 |      1.0 |
| 1459 | 3545 | ATACTGAGA  |  8.750 |    21 |  31 | 0.888889 |     1 |      1.0 |
| 1462 | 3456 | ATACTGAG   |  8.500 |    20 | 193 | 0.864198 |     1 |      1.0 |
| 1458 | 4442 | ATACTGAG   |  8.250 |    24 | 175 | 0.864198 |     1 |      1.0 |
| 5864 | 3961 | TACTGAGA   |  8.250 |    22 | 196 | 0.888889 |     1 |      1.0 |
| 5865 | 4520 | TACTGAGA   |  8.250 |    25 |   2 | 0.888889 |     1 |      1.0 |
| 1463 | 5595 | ATACTGAG   |  7.875 |    29 |  73 | 0.864198 |     1 |      1.0 |
| 5863 | 3546 | TACTGAGA   |  7.750 |    21 |  32 | 0.888889 |     1 |      1.0 |
| 5862 | 4443 | TACTGAG    |  7.250 |    24 | 176 | 0.864198 |     1 |      1.0 |
| 1464 | 5174 | ATACTGA    |  7.125 |    27 | 154 | 0.839506 |     1 |      1.0 |
| 5861 | 5430 | TACTGAG    |  6.875 |    28 | 159 | 0.839506 |     1 |      1.0 |
| 1465 | 4232 | ATACTG     |  6.250 |    23 | 216 | 0.814815 |     1 |      1.0 |

In the context of this small example, it can be useful to employ the
*extendKmers* method to clean up the k-mer output:

```python
extTopTable = sarks.extendKmers(topTable)

```

resulting in extTopTable having the value

|    i |    s | kmer       |   khat | block |  wi |     gini | score | windowed |
|------|------|------------|--------|-------|-----|----------|-------|----------|
| 2257 | 3959 | CATACTGAGA | 10.250 |    22 | 194 | 0.888889 |     1 |      1.0 |
| 2258 | 4518 | CATACTGAGA | 10.250 |    25 |   0 | 0.888889 |     1 |      1.0 |
| 2256 | 3544 | CATACTGAGA |  9.625 |    21 |  30 | 0.864198 |     1 |      1.0 |
| 2257 | 3959 | CATACTGAGA |    NaN |    22 | 194 |      NaN |     1 |      NaN |
| 2258 | 4518 | CATACTGAGA |    NaN |    25 |   0 |      NaN |     1 |      NaN |
| 2256 | 3544 | CATACTGAGA |    NaN |    21 |  30 |      NaN |     1 |      NaN |
| 2259 | 3455 | CATACTGAGA |    NaN |    20 | 192 |      NaN |     1 |      NaN |
| 2255 | 4441 | CATACTGAGA |    NaN |    24 | 174 |      NaN |     1 |      NaN |
| 2257 | 3959 | CATACTGAGA |    NaN |    22 | 194 |      NaN |     1 |      NaN |
| 2258 | 4518 | CATACTGAGA |    NaN |    25 |   0 |      NaN |     1 |      NaN |
| 2260 | 5594 | CATACTGAGA |    NaN |    29 |  72 |      NaN |     1 |      NaN |
| 2256 | 3544 | CATACTGAGA |    NaN |    21 |  30 |      NaN |     1 |      NaN |
| 2255 | 4441 | CATACTGAGA |    NaN |    24 | 174 |      NaN |     1 |      NaN |
| 2261 | 5173 | CATACTGAGA |    NaN |    27 | 153 |      NaN |     1 |      NaN |
| 2254 | 5428 | CATACTGAGA |    NaN |    28 | 157 |      NaN |     1 |      NaN |
| 2262 | 4231 | CATACTGAGA |    NaN |    23 | 215 |      NaN |     1 |      NaN |

The extendKmers method tends to be most useful on small data sets.


