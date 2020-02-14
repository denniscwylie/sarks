orientSeqs <- function(seqs, representative, k=4) {
    if (length(seqs) == 1) {return(seqs)}
    names(seqs) <- seqs
    seqKmers <- kmerCounts(k, seqs, overlap=TRUE)
    revComps <- as.character(Biostrings::reverseComplement(
            Biostrings::DNAStringSet(seqs)))
    revKmers <- kmerCounts(k, revComps, overlap=TRUE)
    reverse <- rep(FALSE, length(seqs))
    repKmers <- as.numeric(seqKmers[representative, ])
    for (i in 1:length(seqs)) {
        if (sum(as.numeric(revKmers[i, ]) * repKmers) >
            sum(as.numeric(seqKmers[i, ]) * repKmers)) {
            reverse[i] = TRUE
        }
    }
    return(ifelse(reverse, revComps, seqs))
}
