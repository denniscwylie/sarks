% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sarks.R
\name{regexCounts}
\alias{regexCounts}
\title{Count occurrences of regular expression}
\usage{
regexCounts(regex, seqs, overlap = FALSE)
}
\arguments{
\item{regex}{character vector of regular expressions to search for}

\item{seqs}{character vector or XStringSet of sequences in which to
search for and count occurrences of regex}

\item{overlap}{logical value: should overlapping occurrences of
regex be counted as multiple hits?}
}
\value{
if length(regex) is one, returns integer vector of counts;
    if length(regex) is more than one, returns matrix of counts:
    one row per sequence in seqs, one column per expression in
    regex
}
\description{
Counts how often a regular expression (or vector of regular
expressions) occurs in each element of a character vector.
}
\examples{
data(simulatedSeqs)
reCounts1 <- regexCounts('AAAAA|TTTTT', simulatedSeqs)
reCounts2 <- regexCounts(c('AAAAA|TTTTT', 'CG'), simulatedSeqs)

}
