% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sarks.R
\name{permutationThresholds}
\alias{permutationThresholds}
\title{Set smoothed score thresholds based on permutation distribution}
\usage{
permutationThresholds(filters, permDist, nSigma = 4)
}
\arguments{
\item{filters}{output from sarksFilters function indicating what
combinations of filter parameters halfWindow, spatialLength,
and minGini to use.}

\item{permDist}{output from permutationDistribution function.}

\item{nSigma}{number of standard deviations above mean of
permutation distribution at which to set threshold for either
windowed or spatially-windowed score.}
}
\value{
named list with two elements: `theta' containing a
    data.frame with the threshold information and `.java'
    containing an R representation of the java object with this
    information.
}
\description{
Calculate thresholds for SArKS k-mer calling from permutation distribution.
}
\examples{
data(simulatedSeqs, simulatedScores)
sarks <- Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
filters <- sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
permDist <- permutationDistribution(sarks, 250, filters, seed=123)
thresholds <- permutationThresholds(filters, permDist, nSigma=2.0)

}
\references{
Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
    SArKS: de novo discovery of gene expression regulatory motif
    sites and domains by suffix array kernel smoothing,
    Bioinformatics, Vol. 35(20), 3944-3952

    \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
}
