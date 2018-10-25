#!/usr/bin/env Rscript

library(philentropy)

## -----------------------------------------------------------------
args = commandArgs(TRUE)
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

## -----------------------------------------------------------------
peaks = read.table(peaksFile,
                   sep='\t', header=TRUE, row.names=NULL,
                   stringsAsFactors=FALSE)

kmers = unique(peaks$kmer)
kmers = unique(sapply(kmers, function(km) {min(km, dnaRevComp(km))}))
names(kmers) = kmers

## -----------------------------------------------------------------
tetraCounts = kmerCounts(kmers, 4, overlap=TRUE)
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
tetraCountsUnorient = data.frame(tetraCountsUnorient, row.names=kmers)
tetraCountsUnorient = tetraCountsUnorient[ , colSums(tetraCountsUnorient) > 0]

## -----------------------------------------------------------------
d = suppressMessages(distance(tetraCountsUnorient, method='jaccard'))
rownames(d) = colnames(d) = rownames(tetraCountsUnorient)
diag(d) = Inf
d[d == 0] = min(d[d > 0]) / 2
diag(d) = 0

## -----------------------------------------------------------------
if (is.na(nClusters)) {
    mds1 = cmdscale(d, k=1)
    mdsKmers = mds1[order(mds1[ , 1]), , drop=FALSE]
    mdsKmers = structure(rownames(mdsKmers), names=rownames(mdsKmers))
    writeLines(mdsKmers)
} else {
    suppressMessages(library(msa))
    hcout = hclust(as.dist(d), method='average')
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
}
