Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
multi-motif domains (MMDs)
(https://www.biorxiv.org/content/early/2018/10/25/133934) whose
presence or absence within each of a set of sequences (contained
within a fasta file) is correlated with numeric scores assigned to
those sequences (contained with a two column tab-separated-values
file, or tsv).

Here we take as an example the downstream mo2015 sequences:
- mo2015\_downstream\_seqs.fa
- mo2015_scores.tsv

Let's try a single parameter set first:
use smallest half-width kappa=250 from paper (-w 250),
use lambda=0 for no spatial smoothing (-l 0),
select Gini impurity filter using gamma=1.1 (-g 1.1)
(this may take a few minutes, depending on speed of computer):

```bash
python3 sarkselect.py -f mo2015_downstream_seqs.fa\
                      -s mo2015_scores.tsv\
                      -w 250\
                      -l 0\
                      -g 1.1\
                      -r 25\
                      -z 4\
                      -o mo2015_downstream_w250_l0_selection
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
- **-o**  
  output directory to be created/overwritten

only used -r 25 permutations for setting threshold for speed here.

- **NOTE:**  
  generally recommend using more permutations (100+) for more stable
  threshold estimation, though smaller numbers can be useful as
  'quick-and-dirty' first look.  Because of the less precise threshold
  estimation, lower permutation numbers can sometimes result in larger
  false positive rates when tested with sarkstest.py below

If we just want to inspect the set of k-mers selected, can use
included z-ranking script; takes one argument: output directory from
sarkselect.py:

```bash
python3 zrank_singly_smoothed_kmers.py mo2015_downstream_w250_l0_selection |\
 column -t
```

| kmer     | halfWindow | minGini |               zmax |
|----------|------------|---------|--------------------|
| TGACCTTG |        250 |     1.1 | 15.923038914859138 |
| TGACCTT  |        250 |     1.1 | 13.578179025162164 |
| GACCTTGG |        250 |     1.1 | 11.179823044834727 |
| GACCTTG  |        250 |     1.1 | 10.385833334059273 |
| AAGGTCA  |        250 |     1.1 |  7.276183749374138 |
| TGTCCTTG |        250 |     1.1 |  5.534467147765855 |
| ACCTTGG  |        250 |     1.1 |  5.392361110568699 |
| TGACCT   |        250 |     1.1 |  4.354476190422092 |
| TGTCCTT  |        250 |     1.1 |  4.270827232664953 |
| CAAGGTC  |        250 |     1.1 |  4.173047347595446 |

This shows each unique k-mer for which a peak was found without
spatial smoothing (hence 'singly_smoothed' in the name of the script)
along with the highest value of z multiplying the
standard-deviation-above-mean (of maximum smoothed suffix scores
obtained after randomly permuting scores assigned to sequences) for
which the peak would have been called.

- **NOTE:**  
  zmax values will vary somewhat b/c permutations are stochastic;
  possible that lowest k-mer may fall below 4 and thus not show up

When spatial smoothing is employed, the process for merging adjacent
peaks makes disrupts this method of motif ranking. I recommend always
including un-spatially-smoothed analyses (lambda=0) as a complement
even when spatially smoothed results are also desired, since some
methods for interpretation of results are easier for lambda=0.

For this parameter combination applied to the downstream sequences,
all the k-mers look pretty similar. Can confirm this by aligning them
into 1 cluster using included cluster\_seqs.R script (NOTE:
cluster\_seqs.R requires *philentropy* and *msa* R libraries):

```bash
Rscript cluster_seqs.R mo2015_downstream_w250_l0_selection/peaks.tsv\
                       1
```

- first argument is tsv file containing a kmer column, such as peaks.tsv,
- second argument is number of clusters to divide kmers into;
- **NOTE:** some kmers will be reverse-complemented in order to align them!)

Output of cluster_seqs.R script should be similar to

```
Cluster 1: CAAGGTCA
--AAGGTCA
---AGGTCA
--AAGGACA
-CAAGGACA
-CAAGGTCA
CCAAGGT--
CCAAGGTC-
-CAAGGTC-
```

(The cluster is named after the k-mer whose average similarity to
all of the other k-mers in the cluster is greatest.)

If you want a position frequency matrix (PFM) representing this alignment:

```bash
Rscript cluster_seqs.R --pfm\
                       mo2015_downstream_w250_l0_selection/peaks.tsv\
                       1
```

Cluster 1: CAAGGTCA

| pos | A | C | G | T |
|-----|---|---|---|---|
|   0 | 0 | 2 | 0 | 0 |
|   1 | 0 | 5 | 0 | 0 |
|   2 | 7 | 0 | 0 | 0 |
|   3 | 8 | 0 | 0 | 0 |
|   4 | 0 | 0 | 8 | 0 |
|   5 | 0 | 0 | 8 | 0 |
|   6 | 2 | 0 | 0 | 6 |
|   7 | 0 | 7 | 0 | 0 |
|   8 | 5 | 0 | 0 | 0 |

Alternately, aligning into 2 clusters yields clusters w/o mismatches:

```bash
Rscript cluster_seqs.R mo2015_downstream_w250_l0_selection/peaks.tsv\
                       2
```

```
Cluster 1: CAAGGTCA
--AAGGTCA
---AGGTCA
-CAAGGTCA
-CAAGGTC-
CCAAGGT--
CCAAGGTC-

Cluster 2: AAGGACA
-AAGGACA
CAAGGACA
```

Selection of the "right" number of clusters is not an easy problem;
it may not even be a well-defined problem depending on the sequences
being analyzed!

Whether you prefer to cluster or not, you can count occurrences of
the selected k-mers or clusters using count_kmers.py:

...for kmers: use -k to specify file with kmer column, -f for fasta

```bash
python3 count_kmers.py -k mo2015_downstream_w250_l0_selection/peaks.tsv\
                       -f mo2015_downstream_seqs.fa >\
                       mo2015_downstream_w250_l0_selection/kmer_counts.tsv
```
...for clusters:

- first save clusters to file:  

  ```bash
  Rscript cluster_seqs.R mo2015_downstream_w250_l0_selection/peaks.tsv\
                         2 >\
                         mo2015_downstream_w250_l0_selection/clusters_2.txt
  ```

- ...and now count: use -c to specify clusters file, -f for fasta
  ```bash
  python3 count_kmers.py -c mo2015_downstream_w250_l0_selection/clusters_2.txt\
                         -f mo2015_downstream_seqs.fa >\
                         mo2015_downstream_w250_l0_selection/cluster_counts.tsv
  ```

- **count_kmers.py:**
  - by default counts both forward and reverse strand matches
    (can get forward matches only using -d (directional) flag)
  - disallows overlapping matches by default
    (can count overlapping kmer matches using -o (overlap) flag)
  - for clusters, it counts a match to any kmer in the cluster
    (cluster matching not compatible with -o option)
  - returns matrix with one row per sequence, one column per kmer/cluster
  - *unless* you pass -l (locate) flag, in which case
     you get one line per match indicating
     1. which sequence matched,
     2. where in the sequence the match occurred (0-based indexing), and
     3. what k-mer or cluster was matched.

This will take a bit, but can once again test false positive rate
- will likely be low here because we used stiff value of z=4 with only
  one parameter combination tested,
- though use of only R=25 permutations in selection above will
  increase variation in results.

```bash
python3 sarkstest.py -f mo2015_downstream_seqs.fa\
                     -s mo2015_scores.tsv\
                     -w 250\
                     -l 0\
                     -g 1.1\
                     -r 250\
                     -z 4\
                     -i mo2015_downstream_w250_l0_selection\
                     -o mo2015_downstream_w250_l0_testing
```

Once you've run this once, you can quickly extract results again via

```bash
python3 sarkstest.py -z 4\
                     -i mo2015_downstream_w250_l0_selection\
                     -o mo2015_downstream_w250_l0_testing
```

Resulting estimate of false positive rate should be similar to:

```
1 / 250
point estimate: 0.4%
95% CI: (0.01%, 2.208%)

```
