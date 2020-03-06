#' Simulated sequences from SArKS paper.
#'
#' 30 simulated DNA sequences (each 250 characters in length) used to
#' illustrate suffix array kernel smoothing method in
#' \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}.
#' Last 10 sequences have the 10-mer CATACTGAGA inserted into them at
#' random locations. Sequences are named "0"-"29".
#'
#' @format Named character vector.
#'
#' @source \url{https://github.com/denniscwylie/sarks/tree/master/examples/simulated_seqs.fa}
"simulatedSeqs"

#' Scores associated with simulated sequences from SArKS paper.
#'
#' Scores associated with simulated DNA sequences used to illustrate
#' suffix array kernel smoothing method in
#' \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}.
#' First 20 sequences assigned score of 0.0, last 10 assigned score of
#' 1.0. Scores are named "0"-"29" so they may be aligned to the
#' corresponding sequences in simulatedSeqs character vector.
#'
#' @format Named numeric vector.
#'
#' @source \url{https://github.com/denniscwylie/sarks/tree/master/examples/simulated_scores.tsv}
"simulatedScores"
