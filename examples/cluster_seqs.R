#!/usr/bin/env Rscript

## library(philentropy)
suppressMessages(library(msa))

## -----------------------------------------------------------------
args = commandArgs(TRUE)
if (all(args %in% c('-h', '--help'))) {
    cat('Usage: cluster_seqs.R [options] <tsv with kmer column> <number of clusters>\nOptions:\n  -d      directional mode\n  -k INT  sub-k-mer length to use in clustering [4]\n')
} else {

if ('-k' %in% args) {
    kargindex = which(args == '-k') + 1
    k = as.integer(args[kargindex])
    args = args[c(-(kargindex-1), -kargindex)]
} else {
    k = 4
}
directional = FALSE
if ('-d' %in% args) {
    directional = TRUE
    args = setdiff(args, '-d')
}
dopfm = FALSE
if ('--pfm' %in% args) {
    dopfm = TRUE
    args = setdiff(args, '--pfm')
}
peaksFile = args[1]
nClusters = if (length(args) > 1) {as.integer(args[2])} else {NA}

## -----------------------------------------------------------------
dnaRevComp = function(s) {
    if (length(s) > 1) {return(sapply(s, dnaRevComp))}
    dnaComplements = c(A='T', C='G', G='C', N='N', T='A')
    return(paste(rev(dnaComplements[ strsplit(s, split='')[[1]] ]),
                 collapse = ''))
}

kmerCounts = function(seq, k, overlap=TRUE) {
    if (is.numeric(k)) {
        nucs = c('A', 'C', 'G', 'T')
        nucs = lapply(1:k, function(...) {nucs})
        k = apply(X = do.call(expand.grid, nucs),
                  MARGIN = 1,
                  FUN = function(x) paste(x, collapse=''))
    }
    names(k) = k
    if (overlap) {
        return(data.frame(lapply(k, function(km) {
            nchar(gsub('[^_]+', '',
                       gsub(paste0('(?=', km, ')'), '_',
                            toupper(seq), perl=TRUE)))
        })))
    } else {
        return(data.frame(lapply(k, function(km) {
            nchar(gsub('[^_]+', '', gsub(km, '_', toupper(seq))))
        })))
    }
}

orientSeqs = function(seqs, representative, k=4) {
    if (length(seqs) == 1) {return(seqs)}
    names(seqs) = seqs
    seqKmers = kmerCounts(seqs, k, overlap=TRUE)
    revComps = sapply(seqs, dnaRevComp)
    revKmers = kmerCounts(revComps, k, overlap=TRUE)
    reverse = rep(FALSE, length(seqs))
    repKmers = as.numeric(seqKmers[representative, ])
    for (i in 1:length(seqs)) {
        if (sum(as.numeric(revKmers[i, ]) * repKmers) >
            sum(as.numeric(seqKmers[i, ]) * repKmers)) {
            reverse[i] = TRUE
        }
    }
    return(ifelse(reverse, revComps, seqs))
}
if (directional) {orientSeqs = function(seqs, ...) {seqs}}

## -----------------------------------------------------------------
if (peaksFile != '-') {
    peaks = read.table(peaksFile,
                       sep='\t', header=TRUE, row.names=NULL,
                       stringsAsFactors=FALSE)
} else {
    f = file('stdin', 'r')
    peaks = data.frame(kmer = unlist(strsplit(readLines(f), ',')),
                       stringsAsFactors = FALSE)
    if (peaks[1, 'kmer'] == 'kmer') {
        peaks = peaks[2:nrow(peaks), , drop=FALSE]
    }
}

kmers = unique(peaks$kmer)
if (!directional) {
    kmers = unique(sapply(kmers, function(km) {min(km, dnaRevComp(km))}))
}
names(kmers) = kmers
kmers = kmers[nchar(kmers) >= k]

## -----------------------------------------------------------------
tetraCounts = kmerCounts(kmers, k, overlap=TRUE)
if (!directional) {
    unorientedTetramers = sort(unique(sapply(
        colnames(tetraCounts),
        function(km) {min(km, dnaRevComp(km))}
    )))

    tetraCountsUnorient = list()
    for (col in unorientedTetramers) {
        if (col == dnaRevComp(col)) {
            tetraCountsUnorient[[col]] = tetraCounts[[col]]
        } else {
            tetraCountsUnorient[[col]] = tetraCounts[[col]] +
                                         tetraCounts[[dnaRevComp(col)]]
        }
    }
    tetraCounts = tetraCountsUnorient
}
tetraCounts = data.frame(tetraCounts, row.names=kmers)
tetraCounts = tetraCounts[ , colSums(tetraCounts) > 0]

## -----------------------------------------------------------------
## d = suppressMessages(philentropy::distance(tetraCounts, method='jaccard'))
pq = as.matrix(tetraCounts) %*%t(as.matrix(tetraCounts))
p2 = rowSums(tetraCounts^2)
d = 1 - pq / (outer(p2, p2, `+`) - pq)
rownames(d) = colnames(d) = rownames(tetraCounts)
diag(d) = Inf
d[d == 0] = min(d[d > 0]) / 2
diag(d) = 0
ddist = as.dist(d)

## -----------------------------------------------------------------
hcout = hclust(ddist, method='average')
if (is.na(nClusters)) {
    suppressMessages(library(cluster))
    integralOptimize = function(f, lower, upper, resolution=10) {
        testPoints = unique(round(seq(lower, upper, by=(upper-lower)/(resolution-1))))
        if (testPoints[length(testPoints)] < upper) {
            testPoints = c(testPoints, upper)
        }
        testVals = sapply(testPoints, f)
        testMax = which.max(testVals)
        newLower = testPoints[ifelse(testMax == 1, testMax, testMax-1)]
        newUpper = testPoints[ifelse(testMax == length(testPoints), testMax, testMax+1)]
        if (newUpper <= (newLower+2)) {
            return(testPoints[testMax])
        } else {
            return(integralOptimize(f, newLower, newUpper, resolution))
        }
    }
    ctouts = data.frame(cutree(hcout, k=2:(nrow(d)-1)), check.names=FALSE)
	nClusters = integralOptimize(
        f = function(nc) {
            mean(cluster:::silhouette.default(ctouts[ , as.character(nc)],
                                              ddist)[ , 3])
        },
        lower = 2,
        upper = nrow(d)-1
    )
}
ctout = cutree(hcout, k=nClusters)
clusters = lapply(1:max(ctout), function(cl) {d[ctout==cl, ctout==cl, drop=FALSE]})
clustReps = sapply(clusters, function(d) {
    rownames(d)[which.min(rowMeans(d))]
})
names(clusters) = clustReps
clustNum = 1
for (clust in clustReps) {
    if (nrow(clusters[[clust]]) > 1) {
        sink('/dev/null')
        msaOut = msa(
            orientSeqs(rownames(clusters[[clust]]), clust),
            type = 'dna'
        )
        sink()
        cat('Cluster ', clustNum, ': ', clust, '\n', sep='')
        if (dopfm) {
            pfm = data.frame(t(consensusMatrix(msaOut@unmasked)[1:4, ]))
            pfm = data.frame(pos=0:(nrow(pfm)-1), pfm)
            write.table(pfm, sep='\t', quote=FALSE, row.names=FALSE)
        } else {
            writeLines(as.character(msaOut@unmasked))
        }
        cat('\n')
    } else {
        if (dopfm) {
            cat('Cluster ', clustNum, ': ', clust, '\n', sep='')
            pfm = data.frame(matrix(0, nrow=nchar(clust), ncol=5))
            colnames(pfm) = c('pos', 'A', 'C', 'G', 'T')
            pfm$pos = 0:(nrow(pfm)-1)
            for (i in 1:nchar(clust)) {
                pfm[i, substr(clust, i, i)] = 1
            }
            write.table(pfm, sep='\t', quote=FALSE, row.names=FALSE)
            cat('\n')
        } else {
            cat('Cluster ', clustNum, ': ', clust, '\n', clust, '\n\n', sep='')
        }
    }
    clustNum = clustNum + 1
}
if (any(ctout == 0)) {
    for (kmi in which(ctout == 0)) {
        clust = rownames(d)[kmi]
        cat('Cluster ', clustNum, ': ', clust, '\n', clust, '\n\n', sep='')
        clustNum = clustNum + 1
    }
}
}
