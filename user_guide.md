# sarks user-guide

Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
multi-motif domains (MMDs)
(https://www.biorxiv.org/content/early/2018/10/25/133934) whose
presence or absence within each of a set of sequences (contained
within a fasta file) is correlated with numeric scores assigned to
those sequences (contained with a two column tab-separated-values
file, or tsv).

For installation instructions, consult [INSTALL.md](INSTALL.md).

# sarks R package

The sarks vignette is the best place to start to learn how to use the
R version of sarks.

The full vignette is available as a pdf if you use the
`"build_vignettes=TRUE"` option when installing sarks in R; otherwise,
you can take a look at the [abridged markdown vignette](sarks_vignette.md).

# sarks command-line execution

Included in the sarks project examples/ folder are several scripts for
usage at the command line:

- [java -jar sarks.jar **select**](#java--jar-sarksjar-select)__
  main step for motif discovery
- [java -jar sarks.jar **test**](#java--jar-sarksjar-test)__
  estimate false positive rate associated with motif set determined
  by **select** command above
- [**extract_kmers.py**](#extract\_kmerspy)  
  extract list of unique k-mers identified by sarks.jar select command
- [**zrank_singly_smoothed_kmers.py**](#zrank\_singly\_smoothed\_kmerspy)  
  rank those k-mers which were selected without use of spatial
  smoothing according to maximum z threshold at which they could be
  identified
- [**cluster_seqs.R**](#cluster\_seqsr)  
  assign k-mers to clusters
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
- [examples/mo2015_upstream_example.md](examples/mo2015_upstream_example.md)  
  this example uses longer sequences and includes spatial smoothing,
  hence it takes longer and requires more RAM.

## java -jar sarks.jar select

Main input for **select** step is provided by the first two arguments
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

- **-e**__
  fixes random number generator seed to generate reproducible output  
  default: no seed
- **-t**__
  number of threads to use during permutational analyses

## java -jar sarks.jar test

The **test** step of sarks requires output from **select** above.

Main arguments -f (fasta) and -s (scores) for **test** step are
identical to those of **select** step above and should be passed same
values. Additional required arguments:

- **-r**  
  number R_2 of permutations used for false positive rate testing,
  need not be same as number R of permutations used to set threshold  
  default: 100
- **-i**  
  input directory for **test** is output directory from **select**;
  needed to determine selection thresholds for testing.
  
Additional optional flags/arguments:

- **-e**__
  fixes random number generator seed to generate reproducible output,
  should take a different value from the seed used in **select** step.__
  default: no seed
- **-t**__
  number of threads to use during permutational analyses

Aside from the output saved to the specified output directory,
**test** command also emits a single number to stdout:
- how many permutations led to non-empty motif sets using the
  thresholds determined by the **select** step.
  - in conjunction with total number R_2 of permutations run during
    **test** step, this can be used to estimate a binomial confidence
    interval for false positive rate.

## extract_kmers.py

Takes one positional argument: tab-separated values (tsv) file with a
header row in which one column is named "kmer" (no quotation marks
required in the column name). The files peaks.tsv or merged_peaks.tsv
(if spatial smoothing was employed) produced by the sarks.jar select
command (in the resulting output directory) are good candidates.

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
directory from **select** step.

One approach to ranking the set of k-mers identified by SArKS is
according to what the maximum value of the -z threshold parameter in
**select** step at which they would be identified.

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

**NOTE:** zmax values from repeated runs of **select** step will vary
somewhat (unless random number generator seed is fixed) because
permutations used to assess null distribution are stochastic.

## cluster_seqs.R

**NOTE:** requires R installation including the libraries:
- **msa** (https://bioconductor.org/packages/release/bioc/html/msa.html)
- **cluster** (https://cran.r-project.org/web/packages/cluster/index.html)

Takes two positional arguments:

1. tab-separated values (tsv) file with a header row in which one
   column is named "kmer" (no quotation marks required in the column
   name), such as the files peaks.tsv or merged_peaks.tsv (if spatial
   smoothing was employed) produced by **select** step (in the
   resulting output directory).
2. number of clusters desired (optional).

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
