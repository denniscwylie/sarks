# Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
# multi-motif domains (MMDs)
# (https://www.biorxiv.org/content/early/2018/10/25/133934) whose
# presence or absence within each of a set of sequences (contained
# within a fasta file) is correlated with numeric scores assigned to
# those sequences (contained with a two column tab-separated-values
# file, or tsv).

# Here we take as an example the upstream mo2015 sequences:
#  mo2015_upstream_seqs.fa
#  mo2015_scores.tsv

# NOTE: the upstream sequences are three times larger than the
#       downstream sequences; also we will be using spatial smoothing
#       for this example. For these reasons, this SArKS example will
#       take both more run-time (hours) and more memory than the
#       downstream example: I would recommend using a machine with at
#       least 48 gigabytes of RAM available.

# This time let's try a few parameter combinations:
# use half-widths kappa=500,2500 from paper (-w 500,2500),
# use lambda=0,100 for either no spatial smoothing or 10 bases (-l 0,10),
# select Gini impurity filter using gamma=1.1 (-g 1.1)
python3 sarkselect.py -f mo2015_upstream_seqs.fa\
                      -s mo2015_scores.tsv\
                      -w 500,2500\
                      -l 0,10\
                      -g 1.1\
                      -r 50\
                      -z 4\
                      -o mo2015_upstream_w500-2500_l0-10_selection
# The arguments here are:
## -----------------------------------------------------------------
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
## -----------------------------------------------------------------
# only used -r 50 permutations for setting threshold for speed here;
# will still take a while (a couple of hours depending on machine used).
# NOTE: generally recommend using more permutations (100+)
#       for more stable threshold estimation, though smaller numbers
#       can be useful as 'quick-and-dirty' first look.
#       Because of the less precise threshold estimation, lower permutation
#       numbers can sometimes result in larger false positive rates
#       when tested with sarkstest.py below

# use zrank_singly_smoothed_kmers.py to rank un-spatially-smoothed k-mers;
# takes one argument: output directory from sarkselect.py:
python3 zrank_singly_smoothed_kmers.py mo2015_upstream_w500-2500_l0-10_selection |\
 column -t
## -----------------------------------------------------------------
# kmer      halfWindow  minGini  zmax
# CAGCCTGG  2500        1.1      5.494459377263551
# CCTGGAA   2500        1.1      5.294138133161724
# CAGCCTG   2500        1.1      5.216757343625498
# CCTGGAAC  2500        1.1      5.11431098357672
# CTGGAACT  2500        1.1      5.058898492086731
# ACCTGC    2500        1.1      4.4586489509093585
# CACCTGC   2500        1.1      4.440212654765286
# TGGCCTG   2500        1.1      4.397117363293189
# CACCTG    2500        1.1      4.366021342281862
# TCCAGC    2500        1.1      4.225106667938478
# CTGGAAC   2500        1.1      4.202105388841035
# CACCTGG   2500        1.1      4.160753168863357
# CTGGTCTA  2500        1.1      4.123361839483047
# GTCCTG    2500        1.1      4.097998587735451
# CTGGAAG   2500        1.1      4.085612733875691
# TTCCAGC   2500        1.1      4.049845865130561
## -----------------------------------------------------------------
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
# determination of significance threshold.
## -----------------------------------------------------------------
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
## -----------------------------------------------------------------
# (Each cluster is named after the k-mer whose average similarity to
#  all of the other k-mers in the cluster is greatest.)
# 
# Feel free to adjust the number of clusters
# (optional second argument to cluster_seqs.R; I omitted above
#  thereby allowing cluster_seqs.R to estimate optimal value
#  using silhouette method).

# Because of the use of spatial smoothing, the z-ranked k-mers are
# only a small fraction of the total k-mer set detected in this
# example. The eighth column of the file 
#  mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv
# contains the final k-mer sets M_spatial (Eq (S16) from the paper)
# broken out by halfWindow (kappa), spatialLength (lambda),
# and minGini (gamma).
#
# We can extract and count the unique k-mers identified using
# extract_kmers.py:
python3 extract_kmers.py\
 mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv |\
 wc -l
# I got 865, though this will vary somewhat from run to run.
# By default, extract_kmers.py considers reverse complements to be
# equivalent (and prints only the lexicographically lower value
# of the identified k-mer and its reverse complement).
# By passing the -d (directional) flag to extract_kmers.py, this
# behavior can be surpressed.

# Can cluster this whole merged k-mer set:
Rscript cluster_seqs.R\
 mo2015_upstream_w500-2500_l0-10_selection/merged_peaks.tsv\
 100 >\
 mo2015_upstream_w500-2500_l0-10_selection/merged_kmer_clusters_100.txt
# this may take a minute.
#
# You can take a look at the resulting merged_kmer_clusters_100.txt
# file to get a better sense of the output; I would suggest starting
# by finding the top few z-ranked k-mers---keep in mind they may be
# reverse-complemented!---among the clusters. For instance, I obtained
# the cluster
## -----------------------------------------------------------------
# Cluster 20: CAGGTGG
# --CAGGTGGT
# --CAGGTGG-
# --CAGGTG--
# -CCAGGTG--
# -GCAGGTG--
# GGCAGGTG--
# --CAGGTGC-
# -GCAGGTGC-
# -GCAGGTGG-
# --GAGGTGG-
# --GAGGTGGA
# ---AGGTGGA
# --CAGGTGGA
## -----------------------------------------------------------------
# containing the reverse complement GCAGGTGG of CCACCTGC.
#
# I chose the cluster size of 100 here to keep the average cluster
# size on the order of 10 for ease of viewing, as we currently use
# clustering mainly to help organize SArKS k-mer results for more
# efficient inspection. It may be useful to try a few different values
# for the number of clusters (or rely on average silhouette estimation
# by omitting number of clusters argument).

# This will take a few hours, but can once again test false positive rate:
python3 sarkstest.py -f mo2015_upstream_seqs.fa\
                     -s mo2015_scores.tsv\
                     -w 500,2500\
                     -l 0,10\
                     -g 1.1\
                     -r 250\
                     -z 4\
                     -i mo2015_upstream_w500-2500_l0-10_selection\
                     -o mo2015_upstream_w500-2500_l0-10_testing
# After you've run this once, you can quickly extract results again via
python3 sarkstest.py -z 4\
                     -i mo2015_upstream_w500-2500_l0-10_selection\
                     -o mo2015_upstream_w500-2500_l0-10_testing
# Resulting estimate of false positive rate should be similar to:
## -----------------------------------------------------------------
# 5 / 250
# point estimate: 2.0%
# 95% CI: (0.653%, 4.605%)
## -----------------------------------------------------------------
# Again, the exact numbers obtained may vary a bit. Increasing the
# number of permutations used in sarkselect.py above will result in
# more tightly controlled false positive rates,
# while increasing the number of permutations used in sarkstest.py will 
# result in more precise (tighter CI) estimates of false positive rate.
