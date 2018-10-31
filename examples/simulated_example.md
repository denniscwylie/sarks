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

Now let's see if SArKS can find this motif using sarkselect.py,
dumping output (-o) to a new directory named sarks\_simulated\_selection

```bash
python3 sarkselect.py -f simulated_seqs.fa\
                      -s simulated_scores.tsv\
                      -w 4\
                      -l 0\
                      -g 1.1\
                      -p\
                      -r 100\
                      -z 2\
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
- **-p**  
  takes no argument; a logical flag indicating to remove redundant kmer output,
  as described in supplementary Section S2.3 of the paper.
  Most useful with small data sets like this example.
- **-r**  
  number R of permutations to use in setting significance thresholds
  for peak-calling
- **-z**  
  multiple z of standard deviations above mean (of maximum smoothed suffix
  scores obtained after randomly permuting scores assigned to sequences)
  defining threshold as described in
  supplementary Section S2.6, Eq (S24-S25) of paper
- **-o**  
  output directory to be created/overwritten

once done, look at what's in new directory:

```bash
ls sarks_simulated_selection
```

should see:

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
be used by the sarkstest.py script below.

Let's look at the peaks.tsv file (cutting off last few columns):

```bash
cut -f 1-9 sarks_simulated_selection/peaks.tsv |\
 sed 's/\t\t/\t-\t/g' |\
 column -t
```

|    i |    s | kmer       | khat  | block |  wi | gini               | score | windowed |
|------|------|------------|-------|-------|-----|--------------------|-------|----------|
| 2255 | 4441 | CATACTGAGA | -     |    24 | 174 | -                  |     1 |          |
| 2259 | 3455 | CATACTGAGA | -     |    20 | 192 | -                  |     1 |          |
| 2260 | 5594 | CATACTGAGA | -     |    29 |  72 | -                  |     1 |          |
| 2256 | 3544 | CATACTGAGA | 9.625 |    21 |  30 | 0.8641975309000001 |     1 |      1.0 |
| 2257 | 3959 | CATACTGAGA | 10.25 |    22 | 194 | 0.8888888889       |     1 |      1.0 |
| 2258 | 4518 | CATACTGAGA | 10.25 |    25 |   0 | 0.8888888889       |     1 |      1.0 |

The khat and gini values for the first three peaks shown are
missing because these peaks have been modified by the extendKmers
method (supplementary Section S2.3.2 of paper).
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
python3 sarkstest.py -f simulated_seqs.fa\
                     -s simulated_scores.tsv\
                     -w 4\
                     -l 0\
                     -g 1.1\
                     -p\
                     -r 1000\
                     -z 2\
                     -i sarks_simulated_selection\
                     -o sarks_simulated_testing
```

The arguments above are identical as for the original selection
step, except for:

- **-r**  
  number R_2 of permutations used for testing,
  need not be same as number R of permutations used to set threshold
- **-i**  
  input directory for testing is output directory from selection step;
  needed to determine selection thresholds for testing
- **-o**  
  output directory for testing; will be created/overwritten

The script sarkstest.py also emits a bit of text on stdout,
should be something like:

```
27 / 1000
point estimate: 2.7%
95% CI: (1.787%, 3.904%)
```

This tells us that 27 out of 1000 permutations generated nonempty
motif sets using the thresholds we calculated using sarkselect.py,
yielding an estimated false positive rate of 2.7% with the
indicated confidence interval, calculated using the Clopper-Pearson
method:

https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval

As the permutation testing can take a while (especially for more
realistically sized data sets), can also re-run sarkstest.py just
to get these false positive rate estimates without recomputing
permutations once they've been generated:

```bash
python3 sarkstest.py -z 2\
                     -i sarks_simulated_selection\
                     -o sarks_simulated_testing
```

which should give same estimates as you obtained above.
