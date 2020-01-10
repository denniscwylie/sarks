test_Sarks <- function() {
    sarks <- Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
    yhat <- sarks$getYhat()
    checkEquals(length(yhat), 7530)
    checkEqualsNumeric(median(yhat), 1/3, tolerance=1e-6)
    checkEqualsNumeric(max(yhat), 1, tolerance=1e-6)
}

test_permutationDistribution <- function() {
    sarks <- Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
    filters <- sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
    permDist <- permutationDistribution(sarks, 20, filters, seed=123)
    checkTrue(is.data.frame(permDist$windowed))
    checkTrue(!any(is.na(permDist$windowed$max)))
    checkTrue(max(permDist$windowed$max) <= 1 + 1e-6)
}

test_kmerPeaks <- function() {
    sarks <- Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
    filters <- sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
    permDist <- permutationDistribution(sarks, 250, filters, seed=123)
    thresholds <- permutationThresholds(filters, permDist, nSigma=2.0)
    peaks <- kmerPeaks(sarks, filters, thresholds)
    checkTrue(all(grepl('TACTGA', peaks$kmer)))
}

test_mergedKmerSubPeaks <- function() {
    sarks <- Sarks(simulatedSeqs, simulatedScores, 4, 3, 1)
    filters <- sarksFilters(halfWindow=4, spatialLength=3, minGini=1.1)
    permDist <- permutationDistribution(sarks, 250, filters, seed=123)
    thresholds <- permutationThresholds(filters, permDist, nSigma=4.0)
    mergedSubPeaks <- mergedKmerSubPeaks(sarks, filters, thresholds)
    checkTrue(all(grepl('CTGAG', mergedSubPeaks$kmer)))
}

test_clusterKmers <- function() {
    kmers <- c(
        'CAGCCTGG', 'CCTGGAA', 'CAGCCTG', 'CCTGGAAC', 'CTGGAACT',
        'ACCTGC', 'CACCTGC', 'TGGCCTG', 'CACCTG', 'TCCAGC',
        'CTGGAAC', 'CACCTGG', 'CTGGTCTA', 'GTCCTG', 'CTGGAAG', 'TTCCAGC'
    )
    kmcl <- clusterKmers(kmers)
    clLens <- vapply(kmcl, length, 0)
    checkEquals(length(kmcl), 6)
    checkEquals(max(clLens), 7)
}

