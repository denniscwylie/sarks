Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
multi-motif domains (MMDs)
(https://www.biorxiv.org/content/early/2018/10/25/133934) whose
presence or absence within each of a set of sequences (contained
within a fasta file) is correlated with numeric scores assigned to
those sequences (contained with a two column tab-separated-values
file, or tsv).

A simple illustrative example is provided here by the pair of files:
- simulated_seqs.fa
- simulated_scores.tsv
Checking the simulated_scores.tsv with

```bash
cat simulated_scores.tsv
```

should display the scores y for each sequence (or block) b:

| b  | y  |
|----|----|
| 0  | 0  |
| 1  | 0  |
|... |... |
| 19 | 0  |
| 20 | 1  |
|... |... |
| 28 | 1  |
| 29 | 1  |

i.e, the first 20 sequences (b from 0 to 19) each have score y=0,
while the last 10 sequences (b from 20 to 29) each have score y=1.
Use the included count_kmers.py script

```bash
python3 count_kmers.py -k CATACTGAGA -f simulated_seqs.fa
```

to confirm that the 10-mer CATACTGAGA correlates perfectly with y

| index | CATACTGAGA |
|-------|------------|
| 0	    | 0          |
| 1	    | 0          |
| 2	    | 0          |
|...    | ...        |
| 19    | 0          |
| 20    | 1          |
|...    | ...        |
| 28    | 1          |
| 29    | 1          |

Now let's see if SArKS can find this motif using select command,
dumping output (-o) to a new directory named sarks\_simulated\_selection

```bash
java -jar sarks.jar\
     select\
     -f simulated_seqs.fa\
     -s simulated_scores.tsv\
     -w 4\
     -l 0\
     -g 1.1\
     -r 100\
     -z 2\
     -e 123\
     -o sarks_simulated_selection
```

The arguments here are:

- **-f**  
  fasta file containing sequences to analyze
- **-s**  
  scores tsv file, col 1: sequence identifiers matching the fasta file,
  col 2: numeric scores
- **-w**  
  half window width for first kernel smoothing pass (kappa in paper),
  can supply multiple values using commas (no spaces)
- **-l**  
  spatial smoothing length (lambda in paper),
  can supply multiple values using commas (no spaces)
- **-g**  
  parameter for calculation of Gini impurity filter (gamma in paper),
  can supply multiple values using commas (no spaces)
- **-r**  
  number R of permutations to use in setting significance thresholds
  for peak-calling
- **-z**  
  multiple z of standard deviations above mean (of maximum smoothed suffix
  scores obtained after randomly permuting scores assigned to sequences)
  defining threshold as described in
  supplementary Section S2.6, Eq (S24-S25) of paper
- **-e**  
  random number generator seed for reproducibility
- **-o**  
  output directory to be created/overwritten

once done, look at what's in new directory:

```bash
ls sarks_simulated_selection
```

should see:

- merged\_peaks.tsv
- peaks.tsv
- permdists\_spatial\_windowed.tsv
- permdists\_theta.tsv
- permdists\_windowed.tsv.

Each of these files is a tab-separated-values (tsv) table. For now, we
consider only two:

- **sarks_simulated_selection/permdists_theta.tsv**  
  displays the thresholds (theta and spatialTheta) for each
  combination of SArKS parameters run (here only one)
- **sarks_simulated_selection/peaks.tsv**  
  contains information similar to Table 1 in SArKS paper for
  all peaks in sequence- and/or spatially-smoothed-score
  notice in particular the column labeled 'kmer'

The other two files (..._windowed.tsv) contain summary statistics for
the distribution of smoothed scores in each null permutation, and will
be used by the sarks.jar test command below.

Let's look at the peaks.tsv file (cutting off last few columns):

```bash
cut -f 1-9 sarks_simulated_selection/peaks.tsv |\
 sed 's/\t\t/\t-\t/g' |\
 column -t
```

|    i |    s | kmer       |  khat | block |  wi |       gini | score | windowed |
|------|------|------------|-------|-------|-----|------------|-------|----------|
| 1458 | 4442 | ATACTGAG   |  8.25 |    24 | 175 | 0.86419755 |   1.0 |      1.0 |
| 1459 | 3545 | ATACTGAGA  |  8.75 |    21 |  31 |  0.8888889 |   1.0 |      1.0 |
| 1460 | 3960 | ATACTGAGA  |  9.25 |    22 | 195 |  0.8888889 |   1.0 |      1.0 |
| 1461 | 4519 | ATACTGAGA  |  9.25 |    25 |   1 |  0.8888889 |   1.0 |      1.0 |
| 1462 | 3456 | ATACTGAGA  |   8.5 |    20 | 193 | 0.86419755 |   1.0 |      1.0 |
| 1463 | 5595 | ATACTGAG   | 7.875 |    29 |  73 | 0.86419755 |   1.0 |      1.0 |
| 2256 | 3544 | CATACTGAGA | 9.625 |    21 |  30 | 0.86419755 |   1.0 |      1.0 |
| 2257 | 3959 | CATACTGAGA | 10.25 |    22 | 194 |  0.8888889 |   1.0 |      1.0 |
| 2258 | 4518 | CATACTGAGA | 10.25 |    25 |   0 |  0.8888889 |   1.0 |      1.0 |
| 5862 | 4443 | TACTGAG    |  7.25 |    24 | 176 | 0.86419755 |   1.0 |      1.0 |
| 5863 | 3546 | TACTGAGA   |  7.75 |    21 |  32 |  0.8888889 |   1.0 |      1.0 |
| 5864 | 3961 | TACTGAGA   |  8.25 |    22 | 196 |  0.8888889 |   1.0 |      1.0 |
| 5865 | 4520 | TACTGAGA   |  8.25 |    25 |   2 |  0.8888889 |   1.0 |      1.0 |

Some of the columns truncated by cut above include information on
which combination of parameters led to each peak; since we
used only one combination for this demonstration
(...-w 4 -l 0 -g 1.1..., or kappa=4, lambda=0, and gamma=1.1)
there are not very interesting in this case.

We know we used z=2 to determine our peak-calling threshold,
but we may still want an estimate of the false-positive rate
associated with this motif set.

We can calculate this by generating an independent set of
permutations and comparing how often the highest observed
smoothed scores exceed the thresholds we determined with
our original permutation set:

```bash
java -jar sarks.jar\
     test\
     -f simulated_seqs.fa\
     -s simulated_scores.tsv\
     -r 1000\
     -i sarks_simulated_selection\
     -e 321
```
The -f (fasta) and -s (scores) arguments above are identical as for
the original selection step. Additional arguments:

- **-r**  
  number R_2 of permutations used for testing,
  need not be same as number R of permutations used to set threshold
- **-i**  
  input directory for testing is output directory from selection step;
  needed to determine selection thresholds for testing
- **-e**  
  optional random number generator seed for reproducibility;
  *should NOT be same number as used for seed in select step!*

The sarks.jar test also emits a number to stdout, such as:

```
16
```

This tells us that 16 / 1000 permutations generated nonempty motif
sets using the thresholds we calculated using sarks.jar select
command, yielding an estimated false positive rate (FPR) of 1.6%.

Recommend calculating a 95% confidence interval (CI) for the FPR using
a method such as:

https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval

Many implementations of this method are available, including
https://epitools.ausvet.com.au/ciproportion

using which I obtained a 95% CI of (0.92%, 2.59%) using
- Sample size: 1000 (number R_2 of permutations used in test step)
- Number positive: 16 (output from java -jar sarks.jar test ...)
- Confidence interval: 0.95
- CI method: Clopper-Pearson exact
