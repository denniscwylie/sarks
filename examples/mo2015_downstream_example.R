#!/usr/bin/env Rscript

# Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
# multi-motif domains (MMDs)
# (https://www.biorxiv.org/content/early/2018/10/25/133934) whose
# presence or absence within each of a set of sequences (contained
# within a fasta file) is correlated with numeric scores assigned to
# those sequences (contained with a two column tab-separated-values
# file, or tsv).

# Here we take as an example the downstream mo2015 sequences:
# - mo2015_downstream_seqs.fa.gz
# - mo2015_scores.tsv

## Because sarks is implemented using rJava, and because the
## default setting for the java virtual machine (JVM) heap space with
## rJava is quite low, you should initialize java with increased heap
## size *before loading sarks*
options(java.parameters = '-Xmx8G')
library(rJava)
.jinit()

library(sarks)
nThreads <- 4  ## set based on how many processors you have available
               ## and would like to devote to running this script

## initialize sarks object using Sarks constructor:
sarks <- Sarks(
    fasta = 'mo2015_downstream_seqs.fa.gz',
    ## fasta file containing sequences to analyze
    scores = 'mo2015_scores.tsv',
    ## scores tsv file,
    ## - col 1: sequence identifiers matching the fasta file,
    ## - col 2: numeric scores
    halfWindow = 250,
    ## half window width for first kernel smoothing pass (kappa in paper)
    spatialLength = 0,
    ## spatial smoothing length (lambda in paper)
    nThreads = nThreads
)

## specify what range of SArKS parameters to investigate:
filters <- sarksFilters(
    halfWindow = 250,
    ## can specify vector w/multiple halfWindow values to investigate
    ## here; for now stick with same as above (but see
    ## mo2015_upstream_example.R as well)
    spatialLength = 0,
    ## again could specify multiple values for spatialLength
    minGini = 1.1
    ## minGini 1.1 specifies how to calculate Gini impurity filter
    ## value (gamma in paper); see
    ## - Eq (13) of vignette for details
)

## estimate null distribution of maximum smoothed scores when no
## association between sequence and score by permutation:
permDist <- permutationDistribution(
    sarks,
    reps = 100,
    ## sample 100 random permutations to estimate null distribution
    filters = filters,
    ## use halfWindow, spatialLength, and minGini values in filters,
    seed = 123
    ## optional: set random generator seed to obtain reproducible results
)

## set thresholds for k-mer peak calling based on permDist:
thresholds <- permutationThresholds(
    filters = filters,
    permDist = permDist,
    nSigma = 4.0
    ## multiple z of standard deviations above mean (of maximum
    ## smoothed suffix scores obtained after randomly permuting scores
    ## assigned to sequences) defining threshold as described in
    ## - Eq (14-15) of vignette
)

## take a look at thresholds:
thresholds$theta
##   halfWindow minGini spatialLength minSpatialGini     theta spatialTheta
## 1        250     1.1             0            1.1 0.7376174           NA
## -----------------------------------------------------------------------------
## - theta column will be defined for rows with spatialLength == 0
## - spatialTheta column will be defined for rows with spatialTheta > 1

## call k-mer peaks using thresholds:
peaks <- kmerPeaks(sarks, filters=filters, thresholds=thresholds)

## =============================================================================
## when not employing spatial smoothing---which introduces some
## additional subtleties---a convenient way to inspect the k-mer
## output is to rank them by a simple z-score metric:
## - first, calculate estimated mean and standard deviation of
##   maximum kernel-smoothed scores from each permutation:
library(dplyr)
permDistStats <- permDist$windowed %>%
    group_by(halfWindow, spatialLength, minGini, minSpatialGini) %>%
    summarize(mean=mean(max), sd=sd(max))
##   halfWindow spatialLength minGini minSpatialGini  mean     sd
## 1        250             0     1.1            1.1 0.620 0.0294
## -----------------------------------------------------------------------------
## now calculate maximum z-scores associated with each unique identified k-mer:
peaks <- left_join(
    peaks,
    permDistStats,
    by = c('halfWindow', 'minGini', 'spatialLength', 'minSpatialGini')
)
peaks$z <- (peaks$windowed - peaks$mean) / peaks$sd
peaks %>%
    group_by(kmer, halfWindow, minGini) %>%
    summarize(zmax = max(z)) %>%
    arrange(-zmax)
##    kmer     halfWindow minGini  zmax
##    <chr>         <int>   <dbl> <dbl>
##  1 TGACCTTG        250     1.1 17.9 
##  2 TGACCTT         250     1.1 15.3 
##  3 GACCTTGG        250     1.1 12.6 
##  4 GACCTTG         250     1.1 11.8 
##  5 AAGGTCA         250     1.1  8.31
##  6 TGTCCTTG        250     1.1  6.37
##  7 ACCTTGG         250     1.1  6.21
##  8 TGACCT          250     1.1  5.06
##  9 TGTCCTT         250     1.1  4.97
## 10 CAAGGTC         250     1.1  4.86
## 11 GACCTTT         250     1.1  4.42
## -----------------------------------------------------------------------------
## note that all of the zmax values should be at least as large as the
## nSigma parameter specified in the permutationThresholds call above
## (defined as 4.0 in this example).

## =============================================================================
## can cluster output k-mers using clusterKmers:
kmClust <- clusterKmers(
    peaks$kmer,
    directional = FALSE
    ## treat reverse-complement sequences as equivalent for clustering
)
kmClust
## $TGACCTTG
## [1] "AAGGTCA"  "ACCTTGG"  "CAAGGTC"  "GACCTTG"  "GACCTTGG" "GACCTTT"  "TGACCTT" 
## [8] "TGACCTTG" "TGACCT"  
##
## $TGTCCTT
## [1] "TGTCCTT"  "TGTCCTTG"

## -----------------------------------------------------------------------------
## can be useful to align k-mers within a cluster to visualize cluster results
## - when clustering has been done with directional=FALSE,
##   need to pick consistent orientation (forward or reverse)
##   for each k-mer.
source('orientSeqs.R')  ## defines function for picking consistent orientations
## Use orientSeqs function just defined to get consistent orientations
## of all k-mers within each cluster:
kmClustOriented <- lapply(names(kmClust), function(representative) {
    unique(orientSeqs(kmClust[[representative]], representative))
    ## unique call here will shrink clusters containing k-mers
    ## which are reverse-complements of each other!
})
names(kmClustOriented) <- names(kmClust)
## Now use msa bioconductor package to obtain multiple sequence
## alignments of k-mers within each cluster:
library(msa)
print(msa(kmClustOriented[[1]], type='dna'))
## [1] TGACCTT--
## [2] TGACCT---
## [3] TGACCTTG-
## [4] -GACCTTT-
## [5] --ACCTTGG
## [6] -GACCTTGG
## [7] -GACCTTG-
## -----------------------------------------------------------------------------
print(msa(kmClustOriented[[2]], type='dna'))
## [1] TGTCCTT-
## [2] TGTCCTTG

## =============================================================================
## finally, let's obtain an estimate of the false positive rate associated
## with the thresholds used to obtain these results:
fpr <- estimateFalsePositiveRate(
    sarks,
    reps = 100,
    ## here use 100 random permutations to estimate false positive rate;
    ## use larger number to shrink width of confidence interval below
    filters = filters,
    thresholds = thresholds,
    seed = 321
    ## use different seed here than was used in defining permDist above!
    ## (permutation set used for testing thresholds should be independent
    ##  of that used to set thresholds.)
)
## extract 95% confidence interval for false positive rate (FPR)
## from fpr object:
fpr$ci
##   method x   n mean lower      upper
## 1  exact 0 100    0     0 0.03621669
## -----------------------------------------------------------------------------
## this indicates that none of 100 permutations yielded a false positive,
## for a point estimate of 0% FPR, with 95% confidence interval of
## (0%, 3.6%)
