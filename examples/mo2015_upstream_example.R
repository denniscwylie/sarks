#!/usr/bin/env Rscript

# Suffix Array Kernel Smoothing, or SArKS, aims to find motifs and
# multi-motif domains (MMDs)
# (https://www.biorxiv.org/content/early/2018/10/25/133934) whose
# presence or absence within each of a set of sequences (contained
# within a fasta file) is correlated with numeric scores assigned to
# those sequences (contained with a two column tab-separated-values
# file, or tsv).

# Here we take as an example the upstream mo2015 sequences:
# - mo2015_upstream_seqs.fa.gz
# - mo2015_scores.tsv

## Because sarks is implemented using rJava, and because the
## default setting for the java virtual machine (JVM) heap space with
## rJava is quite low, you should initialize java with increased heap
## size *before loading sarks*
options(java.parameters = '-Xmx16G')
library(rJava)
.jinit()

library(sarks)
nThreads <- 4  ## set based on how many processors you have available
               ## and would like to devote to running this script

## initialize sarks object using Sarks constructor:
sarks <- Sarks(
    fasta = 'mo2015_upstream_seqs.fa.gz',
    scores = 'mo2015_scores.tsv',
    ## scores tsv file,
    ## - col 1: sequence identifiers matching the fasta file,
    ## - col 2: numeric scores
    halfWindow = 500,  ## kernel half window width (kappa in paper)
    spatialLength = 0, ## spatial smoothing length (lambda in paper)
    nThreads = nThreads
)

## specify what range of SArKS parameters to investigate:
filters <- sarksFilters(
    halfWindow = c(500, 1000),
    ## here specify vector w/multiple halfWindow values to investigate
    spatialLength = c(0, 10),
    minGini = 1.1
    ## minGini 1.1 specifies how to calculate Gini impurity filter
    ## value (gamma in paper); see Eq (13) of vignette for details
)

## estimate null distribution of maximum smoothed scores when no
## association between sequence and score by permutation:
permDist <- permutationDistribution(
    sarks,
    reps = 50,  ## sample 50 permutations to estimate null distribution
    filters = filters,  ## halfWindow, spatialLength, minGini values
    seed = 123  ## optional: set seed to obtain reproducible results
)

## set thresholds for k-mer peak calling based on permDist:
thresholds <- permutationThresholds(
    filters = filters,
    permDist = permDist,
    nSigma = 4.0
    ## multiple z of standard deviations above mean (of maximum
    ## smoothed suffix scores *or* spatially smoothed scores,
    ## depending on whether spatial smoothing employed) defining
    ## threshold as described in Eq (14-15) of vignette
)

## take a look at thresholds:
thresholds$theta
##   halfWindow minGini spatialLength minSpatialGini     theta spatialTheta
## 1        500     1.1             0            1.1 0.5819998           NA
## 2        500     1.1            10            1.1        NA    0.3689503
## 3       1000     1.1             0            1.1 0.4359220           NA
## 4       1000     1.1            10            1.1        NA    0.2634068
## -----------------------------------------------------------------------------
## - theta column will be defined for rows with spatialLength == 0
## - spatialTheta column will be defined for rows with spatialTheta > 1

## call k-mer peaks using thresholds:
## use mergedKmerSubPeaks function when spatial smoothing is included
## (see discussion in vignette section 1.5)
peaks <- mergedKmerSubPeaks(sarks, filters=filters, thresholds=thresholds)

## =============================================================================
## use clusterKmers to get a better sense of output:
kmClust <- clusterKmers(
    peaks$kmer,
    directional = FALSE
    ## treat reverse-complement sequences as equivalent for clustering
)
## how many clusters were identified?
length(kmClust)   ## 104
## take a look at first cluster:
kmClust[1]
## $CCACCTGC
## [1] "CCACCTGC"  "CCACCTGCC" "CACCTGC"   "CACCTGCC"  "GCAGGTGC"  "CCACCTGCT"
## [7] "TCCACCTGC"

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
## alignments of k-mers within a particular cluster:
library(msa)
print(msa(kmClustOriented[[1]], type='dna'))
## [1] --CACCTGC-
## [2] -GCACCTGC-
## [3] -CCACCTGCC
## [4] --CACCTGCC
## [5] -CCACCTGC-
## [6] -CCACCTGCT
## [7] TCCACCTGC-


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
##   method x   n mean       lower      upper
## 1  exact 1 100 0.01 0.000253146 0.05445939
## -----------------------------------------------------------------------------
## this indicates that 1 out of 100 permutations yielded a false positive,
## for a point estimate of 0.01 (1%) FPR, with 95% confidence interval
## of (0.025%, 5.4%)


## =============================================================================
bi <- blockInfo(sarks, 'ENSMUST00000024702', filters, thresholds)

library(ggplot2)
ggo <- ggplot(bi, aes(x=wi+1))  ## +1 because R indexing is 1-based
ggo <- ggo + geom_point(aes(y=windowed), alpha=0.3, size=0.5)
ggo <- ggo + facet_grid(halfWindow ~ spatialLength)
ggo <- ggo + geom_line(aes(y=spatialWindowed), alpha=1, size=0.5)
ggo <- ggo + geom_hline(aes(yintercept=theta), color='red')
ggo <- ggo + geom_hline(aes(yintercept=spatialTheta), color='red')
ggo <- ggo + ylab('yhat') + ggtitle('Paqr4')
ggo <- ggo + theme_bw()
print(ggo)

sum(bi$spatialWindowed > bi$spatialTheta, na.rm=TRUE)  ## 6
summary(bi[(!is.na(bi$spatialTheta)) &
           (bi$spatialWindowed > bi$spatialTheta), 'gini'])
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.9991  0.9992  0.9993  0.9993  0.9993  0.9993
summary(bi$gini)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.9863  0.9988  0.9989  0.9990  0.9993  0.9995 
