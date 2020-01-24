# Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
# multi-motif domains (MMDs)
# (https://www.biorxiv.org/content/early/2018/10/25/133934) whose
# presence or absence within each of a set of sequences (contained
# within a fasta file) is correlated with numeric scores assigned to
# those sequences (contained with a two column tab-separated-values
# file, or tsv).

# Here we take as an example the upstream mo2015 sequences:
#  mo2015_upstream_seqs.fa.gz
#  mo2015_scores.tsv

# NOTE: the upstream sequences are three times larger than the
#       downstream sequences; also we will be using spatial smoothing
#       for this example. For these reasons, this SArKS example will
#       take both more run-time (minutes) and more memory than the
#       downstream example.

# This time let's try a few parameter combinations:
# use half-widths kappa=500,2500 from paper (-w 500,2500),
# use lambda=0,100 for either no spatial smoothing or 10 bases (-l 0,10),
# select Gini impurity filter using gamma=1.1 (-g 1.1)
java -Xmx8G -jar sarks.jar\
     select\
     -f mo2015_upstream_seqs.fa.gz\
     -s mo2015_scores.tsv\
     -w 500,2500\
     -l 0,10\
     -g 1.1\
     -r 50\
     -z 4\
     -o mo2015_upstream_w500-2500_l0-10_selection\
     -t 4\
     -e 123
# The arguments here are:
## -----------------------------------------------------------------------------
# -f <fasta file containing sequences to analyze>
# -s <scores tsv file, col 1: sequence identifiers matching the fasta file,
#                      col 2: numeric scores>
# -w <half window width for first kernel smoothing pass (kappa in paper),
#     can supply multiple values using commas (no spaces)>
# -l <spatial smoothing length (lambda in paper),
#     can supply multiple values using commas (no spaces)>
# -g <parameter for calculation of Gini impurity filter (gamma in paper),
#     can supply multiple values using commas (no spaces)>
# -r <number R of permutations to use in setting significance thresholds
#     for peak-calling>
# -z <multiple z of standard deviations above mean (of maximum smoothed suffix
#     scores obtained after randomly permuting scores assigned to sequences)
#     defining threshold as described in
#     supplementary Section S2.6, Eq (S24-S25) of paper>
# -o <output directory to be created/overwritten>
# -t <number of threads to use for permutational analyses>
# -e <random number generator seed for reproducibility>
## -----------------------------------------------------------------------------
# only used -r 50 permutations for setting threshold for speed here.
# NOTE: generally recommend using more permutations (100+)
#       for more stable threshold estimation, though smaller numbers
#       can be useful as 'quick-and-dirty' first look.
#       Because of the less precise threshold estimation, lower permutation
#       numbers can sometimes result in larger false positive rates
#       when tested with test command below

# use zrank_singly_smoothed_kmers.py to rank un-spatially-smoothed k-mers;
# takes one argument: output directory from sarks.jar select command:
python3 zrank_singly_smoothed_kmers.py mo2015_upstream_w500-2500_l0-10_selection |\
 column -t
## -----------------------------------------------------------------------------
# kmer      halfWindow  minGini  zmax
# CAGCCTGG  2500        1.1      5.441745218958897
# CCTGGAA   2500        1.1      5.252588923599599
# CAGCCTG   2500        1.1      5.179519880370132
# CCTGGAAC  2500        1.1      5.082784566055976
# CTGGAACT  2500        1.1      5.030460072171956
# ACCTGC    2500        1.1      4.463664742078993
# CACCTGC   2500        1.1      4.446257097956867
# TGGCCTG   2500        1.1      4.405563429072066
# CACCTG    2500        1.1      4.376200775586434
# TCCAGC    2500        1.1      4.2431397501378
# CTGGAAC   2500        1.1      4.221419397730943
# CACCTGG   2500        1.1      4.182371837923922
# CTGGTCTA  2500        1.1      4.147064991934326
# GTCCTG    2500        1.1      4.123114348600267
# CTGGAAG   2500        1.1      4.111419336721907
# TTCCAGC   2500        1.1      4.07764566340248
# CCACCTGC  500         1.1      4.023236346571467
## -----------------------------------------------------------------------------
# As in the downstream example, This shows each unique k-mer for which a
# peak was found without spatial smoothing (hence 'singly_smoothed' in
# the name of the script) along with the highest value of z
# multiplying the standard-deviation-above-mean (of maximum smoothed
# suffix scores obtained after randomly permuting scores assigned to
# sequences) for which the peak would have been called.
#
# NOTE: zmax values will vary somewhat b/c permutations are stochastic;
#       possible that lowest k-mer may fall below 4 and thus not show up
#
# (When spatial smoothing is employed, the process for merging adjacent
#  peaks makes disrupts this method of motif ranking. I recommend always
#  including un-spatially-smoothed analyses (lambda=0) as a complement
#  even when spatially smoothed results are also desired, since some
#  methods for interpretation of results are easier for lambda=0.)

# A quick glance through the above list of k-mers suggests some but
# not all of them are quite similar to each other. In order to cluster
# just this set of kmers, let's save the zrank-ed output...
python3 zrank_singly_smoothed_kmers.py\
 mo2015_upstream_w500-2500_l0-10_selection >\
 mo2015_upstream_w500-2500_l0-10_selection/zranked_singly_smoothed.tsv
# ...so that we can use clusters_seqs.R to cluster the k-mers:
# (cluster_seqs.R can be applied to any tsv file with a column named 'kmer')
Rscript cluster_seqs.R\
 mo2015_upstream_w500-2500_l0-10_selection/zranked_singly_smoothed.tsv
# (first argument is tsv file containing a kmer column, such as peaks.tsv,
#  second argument would be number of clusters to divide kmers into,
#  though here we omit this argument and rely on average silhouette
#  method to estimate optimal number of clusters;
#  NOTE: some kmers will be reverse-complemented in order to align them!)
#
# This should result similar output to that shown below; may be
# slightly different due to stochasticity in permutational
# determination of significance threshold if seed is not set or is
# changed.
## -----------------------------------------------------------------------------
# Cluster 1: CAGCCTGG
# CAGCCTGG
# CAGCCTG-
#
# Cluster 2: CTGGAAC
# CCTGGAA--
# CCTGGAAC-
# -CTGGAACT
# -CTGGAAC-
# -CTGGAAG-
# GCTGGA---
# GCTGGAA--
#
# Cluster 3: CACCTGC
# CACCTG-
# CACCTGG
# -ACCTGC
# CACCTGC
#
# Cluster 4: CAGGCCA
# CAGGCCA
#
# Cluster 5: CTGGTCTA
# CTGGTCTA
#
# Cluster 6: CAGGAC
# CAGGAC
## -----------------------------------------------------------------------------
# (Each cluster is named after the k-mer whose average similarity to
#  all of the other k-mers in the cluster is greatest.)
# 
# Feel free to adjust the number of clusters
# (optional second argument to cluster_seqs.R; I omitted above
#  thereby allowing cluster_seqs.R to estimate optimal value
#  using silhouette method).

# Because of the use of spatial smoothing, the z-ranked k-mers are
# only a small fraction of the total k-mer set detected in this
# example. The third column of the file 
#  mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv
# contains the final k-mer sets M_spatial (Eq (25) from vignette)
# broken out by halfWindow (kappa), spatialLength (lambda),
# and minGini (gamma).
#
# We can extract and count the unique k-mers identified using
# extract_kmers.py:
python3 extract_kmers.py\
 mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv |\
 wc -l
# I got 1237, though this will vary somewhat if you use a different
# seed (or do not set the seed parameter).
# 
# By default, extract_kmers.py considers reverse complements to be
# equivalent (and prints only the lexicographically lower value of the
# identified k-mer and its reverse complement). By passing the -d
# (directional) flag to extract_kmers.py, this behavior can be
# suppressed.

# Can cluster this whole merged k-mer set:
Rscript cluster_seqs.R\
 mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv >\
 mo2015_upstream_w500-2500_l0-10_selection/merged_kmer_clusters.txt
# this may take a minute.
#
# You can take a look at the resulting merged_kmer_clusters.txt file
# to get a better sense of the output; I would suggest starting by
# finding the top few z-ranked k-mers---keep in mind they may be
# reverse-complemented!---among the clusters. For instance, I obtained
# the cluster
## -----------------------------------------------------------------------------
# Cluster 1: GCAGGTGGA
# --GCAGGTGGA
# -CGCAGGTGGA
# GTGCAGGTGGA
# -TGCAGGTGGA
# --GCAGGTGG-
# --GCAGGTGGT
## -----------------------------------------------------------------------------
# containing the reverse complement GCAGGTGG of CCACCTGC.
#
# This may take a few minutes, but can once again test false positive rate:
java -jar sarks.jar\
     test\
     -f mo2015_upstream_seqs.fa.gz\
     -s mo2015_scores.tsv\
     -r 250\
     -i mo2015_upstream_w500-2500_l0-10_selection\
     -t 4\
     -e 321
# Resulting estimate of number of false positives should be similar to:
## -----------------------------------------------------------------------------
# 5
## -----------------------------------------------------------------------------
# This gives a point estimate of false positive rate (FPR) of 5 / 250
# or 2%.
# 
# Using the Clopper-Pearson method for estimating a binomial
# confidence interval with sample size 250 and 5 positives, obtain
# (0.65%, 4.6%)
#
# Again, if you change (or don't set) the random number generator seed
# (-e), the exact numbers obtained may vary a bit. Increasing the
# number of permutations used in the select command above will result
# in more tightly controlled FPRs, while increasing the number of
# permutations used in the test command will result in more precise
# (tighter CI) estimates of FPR.
