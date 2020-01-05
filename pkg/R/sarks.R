.onLoad = function(libname, pkgname) {
    rJava::.jpackage(pkgname, lib.loc=libname)
}

#' Suffix Array Kernel Smoothing
#'
#' Sarks class implements suffix array kernel smoothing for de novo
#' correlative motif discovery.
#' 
#' @param fasta specification of fasta file containing sequences to be
#'     analyzed; may also be a named character vector whose elements
#'     are sequences to be analyzed.
#'
#' @param scores specification of scores associated with sequences in
#'     fasta argument; can be provided as two column tab-delimited
#'     file (should have header, first column should provide sequence
#'     names identical to those in fasta argument, second column
#'     should have numeric scores) or may be a named numeric vector.
#'
#' @param halfWindow half-width of smoothing window (integer).
#'
#' @param spatialLength full length of spatial smoothing window
#'     (integer); use 0 to disable spatial smoothing.
#'
#' @param nThreads number of threads to use for computing permutation distributions.
#'
#' @return R representation of java Sarks object.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#'
#' @export
Sarks = function(fasta, scores, halfWindow, spatialLength=0L, nThreads=1L) {
    if (length(fasta) > 1) {
        seqs = character(2*length(fasta))
        seqs[seq(1, length(seqs)-1, by=2)] = paste0('>', names(fasta))
        seqs[seq(2, length(seqs), by=2)] = fasta
        tmpfasta = tempfile(pattern='seq')
        writeLines(seqs, tmpfasta)
        fasta = tmpfasta
    }
    if (length(scores) > 1) {
        tmpscores = tempfile(pattern='scores')
        utils::write.table(
            data.frame(block=names(scores), score=scores, stringsAsFactors=FALSE),
            tmpscores,
            sep = '\t',
            quote = FALSE,
            row.names = FALSE
        )
        scores = tmpscores
    }
    sarks = rJava::.jnew('dcw/sarks/Sarks', fasta, scores,
                         as.integer(halfWindow), as.integer(spatialLength),
                         as.integer(nThreads), TRUE)
    if (exists('tmpfasta')) {suppressWarnings(file.remove(tmpfasta))}
    if (exists('tmpscores')) {suppressWarnings(file.remove(tmpscores))}
    return(sarks)
}

#' Smoothing window and Gini impurity filter settings
#' 
#' Sarks methodology involves testing a range of different filter
#' parameter values; sarksFilters builds set of filters with all
#' combinations of desired halfWindow, spatialLength, and minGini
#' values.
#' 
#' @param halfWindow integer vector of halfWindow values to test.
#'
#' @param spatialLength integer vector of spatialLength values to
#'     test; use a single 0 value to disable spatial smoothing.
#' 
#' @param minGini numeric vector giving minimum Gini impurity value(s)
#'     for suffix position to be analyzed; use a value above 1 to
#'     calculate minimum Gini impurity based on median of observed
#'     Gini impurities so as to constrain variance under permutation
#'     testing to less than minGini multiples of median variance.
#'
#' @return R representation of java object containing specified
#'     combinations of filter parameters for running permutation
#'     tests.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=c(4, 8), spatialLength=c(0, 5), minGini=1.1)
#'
#' @export
sarksFilters = function(halfWindow, spatialLength, minGini=1.1) {
    filters = rJava::.jnew('java/util/ArrayList')
    for (hw in halfWindow) {for (sl in spatialLength) {for (mg in minGini) {
        pars = rJava::.jnew('java/util/HashMap')
        pars$put('halfWindow', rJava::.jnew('java/lang/Integer', as.integer(hw)))
        pars$put('spatialLength', rJava::.jnew('java/lang/Integer', as.integer(sl)))
        pars$put('minGini', rJava::.jnew('java/lang/Double', mg))
        pars$put('minSpatialGini', rJava::.jnew('java/lang/Double', mg))
        filters$add(pars)
    }}}
    return(filters)
}


#' Estimating distribution of maximum smoothed sequence scores under
#' permutation
#'
#' Run permutation test using the specified number of repetitions,
#' keeping track of maximum observed windowed and spatially-windowed
#' smoothed scores for each combination of filter parameters for each
#' permutation.
#'
#' @param sarks Sarks object to test.
#'
#' @param reps integer specifying how many repetitions to test.
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini to use.
#'
#' @param seed optional seed for random number generator (use in case
#'     reproducibility of output is desired).
#'
#' @return named list with three elements: `windowed' containing a
#'     data.frame with the maximum smoothed scores for each
#'     permutation at each combination of filter parameter values,
#'     `spatial' containing a data.frame with the maximum
#'     spatially-smoothed scores for each permutation and each filter
#'     parameter specification, and `.java' containing the R
#'     representation of the java object encoding this information.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#'
#' @export
permutationDistribution = function(sarks, reps, filters, seed=NULL) {
    jreps = rJava::.jnew('java/lang/Integer', as.integer(reps))
    if (length(seed) == 0) {
        obj = sarks$permutationDistribution(jreps, filters)
    } else {
        obj = sarks$permutationDistribution(jreps, filters, as.integer(seed))
    }
    jw = rJava::J('dcw/sarks/Sarks')$printPermDist(filters, obj, NULL, FALSE)
    js = rJava::J('dcw/sarks/Sarks')$printPermDist(filters, obj, NULL, TRUE)
    out = list(
        windowed = utils::read.table(
            textConnection(jw),
            sep='\t', header=TRUE, row.names=NULL, na.strings='null',
            stringsAsFactors=FALSE, check.names=FALSE
        ),
        spatial = utils::read.table(
            textConnection(js),
            sep='\t', header=TRUE, row.names=NULL, na.strings='null',
            stringsAsFactors=FALSE, check.names=FALSE
        ),
        .java = obj
    )
    return(out)
}


#' Set smoothed score thresholds based on permutation distribution
#'
#' Calculate thresholds for SArKS k-mer calling from permutation distribution.
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini to use.
#'
#' @param permDist output from permutationDistribution function.
#'
#' @param nSigma number of standard deviations above mean of
#'     permutation distribution at which to set threshold for either
#'     windowed or spatially-windowed score.
#'
#' @return named list with two elements: `theta' containing a
#'     data.frame with the threshold information and `.java'
#'     containing an R representation of the java object with this
#'     information.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#'
#' @export
permutationThresholds = function(filters, permDist, nSigma) {
    obj = rJava::J('dcw/sarks/SarksUtilities')$thresholdsFromPermutations(
            permDist$.java, nSigma)
    jThresh = rJava::J('dcw/sarks/Sarks')$printThresholds(filters, obj, NULL)
    return(list(
        theta = utils::read.table(
            textConnection(jThresh),
            sep='\t', header=TRUE, row.names=NULL, na.strings='null',
            stringsAsFactors=FALSE, check.names=FALSE
        ),
        .java = obj
    ))
}


#' Call k-mer peaks
#'
#' SArKS identifies sets of short subsequences (k-mers) whose presence
#' as substrings of sequences from the input sequence set tends to be
#' associated with elevated sequence scores. Such k-mers are
#' identified as ``peaks'' where kernel-smoothed scores exceed
#' specified thresholds (generally set by permutation method).
#'
#' @param sarks Sarks object to use for k-mer peak calling.
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini to use.
#'
#' @param thresholds output from permutationThresholds specifying
#'     thresholds for k-mer peak calling.
#'
#' @param peakify logical value specifying whether to restrict output
#'     to only spatial positions at which the smoothed score is at
#'     least as high as either neighboring position or not.
#'
#' @param kMax integer value indicating the maximum k-mer length to be
#'     reported.
#'
#' @return data.frame containing called k-mer peak information.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#' peaks = kmerPeaks(sarks, filters, thresholds)
#'
#' @export
kmerPeaks = function(sarks, filters, thresholds, peakify=TRUE, kMax=12L) {
    kMax = as.integer(kMax)
    eyes = sarks$filter(filters, thresholds$.java,
                        rJava::.jnew('java/lang/Boolean', peakify))
    jOut = sarks$printPeaks(filters, thresholds$.java, eyes, kMax)
    return(utils::read.table(
        textConnection(jOut),
        sep='\t', header=TRUE, row.names=NULL, na.strings='null',
        stringsAsFactors=FALSE, check.names=FALSE
    ))
}


#' Identify and merge k-mer sub-peaks within multi-motif domains
#'
#' When spatials smoothing is employed, SArKS identifies spatial
#' windows containing elevated spatially-averaged sequence-smoothed
#' scores (multi-motif domains, or MMDs). This function finds k-mers
#' within these MMDs whose sequence-smoothed scores are above the
#' threshold used for MMD calling and merges such k-mers when their
#' spatial positions overlap.
#'
#' @param sarks Sarks object to use for k-mer peak calling.
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini to use.
#'
#' @param thresholds output from permutationThresholds specifying
#'     thresholds for k-mer peak calling.
#'
#' @param peakify logical value specifying whether to restrict initial
#'     k-mer peak calling to only spatial positions at which the
#'     smoothed score is at least as high as either neighboring
#'     position (or not).
#'
#' @param kMax integer value indicating the maximum k-mer length for
#'     initial k-mer peak calling.
#'
#' @return data.frame containing called k-mer peak information.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 3, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=3, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=4.0)
#' mergedSubPeaks = mergedKmerSubPeaks(sarks, filters, thresholds)
#'
#' @export
mergedKmerSubPeaks = function(sarks, filters, thresholds,
                              peakify=TRUE, kMax=12L) {
    kMax = as.integer(kMax)
    eyes = sarks$filter(filters, thresholds$.java,
                        rJava::.jnew('java/lang/Boolean', peakify))
    mergedKmerIntervals = sarks$multiMergeKmerIntervals(
            eyes, thresholds$.java, filters, kMax)
    jOut = sarks$printMergedSubPeaks(filters, thresholds$.java,
                                     mergedKmerIntervals, NULL, kMax)
    return(utils::read.table(
        textConnection(jOut),
        sep='\t', header=TRUE, row.names=NULL, na.strings='null',
        stringsAsFactors=FALSE, check.names=FALSE
    ))
}


#' Estimating SArKS false positive rate
#'
#' Run second permutation test using the specified number of
#' repetitions, keeping track of maximum observed windowed and
#' spatially-windowed smoothed scores for each combination of filter
#' parameters for each permutation, and comparing these values to
#' thresholds determined by first round of permutation testing.
#'
#' @param sarks Sarks object to test.
#'
#' @param reps integer specifying how many repetitions to test.
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini to use.
#'
#' @param thresholds output from permutationThresholds specifying
#'     thresholds for k-mer peak calling.
#'
#' @param seed optional seed for random number generator (use in case
#'     reproducibility of output is desired).
#'
#'     NOTE: do not use the same seed passed to initial
#'     permutationDistribution call used to set thresholds.
#'
#' @param conf.level level of confidence to be used in the false
#'     positive rate confidence interval.
#'
#' @return named list with three elements: `permutation' containing the output from permutationDistribution run 
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#' fpr = estimateFalsePositiveRate(sarks, 250, filters, thresholds, seed=123456)
#'
#' @export
estimateFalsePositiveRate = function(sarks, reps,
                                     filters, thresholds,
                                     seed=NULL, conf.level=0.95) {
    permTest = permutationDistribution(sarks, as.integer(reps), filters)
    falsePositives = rJava::J('dcw/sarks/SarksUtilities')$falsePositives(
        permTest$.java,
        thresholds$.java
    )
    return(list(
        fp = falsePositives,
        ci = binom::binom.exact(falsePositives, reps, conf.level),
        permutations = permTest
    ))
}


#' Prune nested k-mer intervals from called k-mer peak set.
#'
#' Every k-mer identified by SArKS is derived as a substring defined
#' by the interval running position i to position i+k-1 of the
#' concatenation of all input sequences. In some cases a j-mer (with j
#' < k) may be separately identified as a peak by SArKS for which the
#' j-mer interval is entirely contained within [i, i+k-1]; this
#' function removes such nested intervals from the reported collection
#' of peaks.
#'
#' @param intervals data.frame containing called k-mer peaks
#'     information (format as output from kmerPeaks function).
#'
#' @param start name of column in intervals data.frame containing
#'     interval start coordinates
#'
#' @param end name of column in interval data.frame containing
#'     interval end coordinates; if no such column present, default
#'     NULL value indicates that end coordinates should be obtained by
#'     adding nchar(intervals$kmer) to the start coordinates to obtain
#'     end coordinates.
#'
#' @return modified data.frame containing called k-mer peaks
#'     information (format as output from kmerPeaks function).
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#' peaks = kmerPeaks(sarks, filters, thresholds)
#' prunedPeaks = pruneIntervals(peaks)
#' @export
pruneIntervals = function(intervals, start='s', end=NULL) {
    starts = intervals[[start]]
    if (length(end) == 0) {
        ends = starts + nchar(intervals$kmer)
    } else {
        ends = intervals[[end]]
    }
    ncl = IRanges::NCList(IRanges::IRanges(starts, ends))
    nests = IRanges::findOverlaps(ncl, type='within',
                                  drop.self=TRUE, drop.redundant=TRUE)
    if (length(nests@from) > 0) {
        intervals = intervals[-nests@from, ]
    }
    return(intervals)
}


#' Extend k-mers in length where possible
#' 
#' Extend k-mers when adding flanking characters from region in
#' input sequence from which they are derived would result in another
#' reported l-mer string (l > k).
#'
#' @param sarks Sarks object used to obtain k-mer peak call set.
#'
#' @param sarksTable data.frame containing called k-mer peaks
#'     information (format as output from kmerPeaks function).
#'
#' @return modified data.frame containing called k-mer peaks
#'     information (format as output from kmerPeaks function).
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#' peaks = kmerPeaks(sarks, filters, thresholds)
#' prunedPeaks = pruneIntervals(peaks)
#' extendedPeaks = extendKmers(sarks, prunedPeaks)
#'
#' @export
extendKmers = function(sarks, sarksTable) {
    maxLen = max(nchar(sarksTable$kmer))
    oldS = rep(-1, nrow(sarksTable))
    while (!all(oldS == sarksTable$s)) {
        oldS = sarksTable$s
        kmerSet = unique(sarksTable$kmer)
        kmerLens = nchar(sarksTable$kmer)
        for (idx in which(kmerLens < maxLen)) {
            doBreak = FALSE            
            s = sarksTable[idx, 's']
            kmerLen = nchar(sarksTable[idx, 'kmer'])
            shiftInt = c(s, s+kmerLen)
            for (ext in (maxLen-kmerLen):0) {
                for (shift in (-ext):0) {
                    shiftInt = c(s+shift, s+kmerLen+ext+shift)
                    shiftKmer = sarks$kmer(as.integer(shiftInt[1]),
                                           as.integer(shiftInt[2]-shiftInt[1]))
                    if (shiftKmer %in% kmerSet) {
                        sarksTable[idx, 's'] = sarksTable[idx, 's'] + shift
                        sarksTable[idx, 'i'] = sarks$s2i(sarksTable[idx, 's'])
                        sarksTable[idx, 'wi'] = sarksTable[idx, 'wi'] + shift
                        sarksTable[idx, 'kmer'] = shiftKmer
                        sarksTable[idx, 'khat'] = NA
                        sarksTable[idx, 'gini'] = NA
                        sarksTable[idx, 'windowed'] = NA
                        sarksTable[idx, 'spatial_windowed'] = NA
                        doBreak = TRUE
                        break
                    }
                }
                if (doBreak) {break}
            }
        }
    }
    return(sarksTable)
}


#' Get sarks smoothed scores for input sequence
#'
#' Returns a data.frame containing the smoothed scores (including
#' spatially smoothed scores, if applicable) as well as other useful
#' sarks parameters for one or more specified input sequences
#' (blocks).
#'
#' @param sarks Sarks object from which information will be derived.
#'
#' @param block character vector of names of sequence(s) for which
#'     results are desired
#'
#' @param filters output from sarksFilters function indicating what
#'     combinations of filter parameters halfWindow, spatialLength,
#'     and minGini were used.
#'
#' @param thresholds output from permutationThresholds specifying
#'     thresholds used for k-mer peak calling.
#'
#' @param kMax integer value indicating the maximum k-mer length to be
#'     reported.
#'
#' @return data.frame in same format as result of kmerPeaks giving
#'     detailed information about every spatial position within
#'     specified sequences/blocks.
#'
#' @references Wylie, D.C., Hofmann, H.A., and Zemelman, B.V. (2019)
#'     SArKS: de novo discovery of gene expression regulatory motif
#'     sites and domains by suffix array kernel smoothing,
#'     Bioinformatics, Vol. 35(20), 3944-3952
#'
#'     \url{https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797}
#'
#' @examples
#' sarks = Sarks(simulatedSeqs, simulatedScores, 4, 0, 1)
#' filters = sarksFilters(halfWindow=4, spatialLength=0, minGini=1.1)
#' permDist = permutationDistribution(sarks, 250, filters, seed=123)
#' thresholds = permutationThresholds(filters, permDist, nSigma=2.0)
#' bi24 = blockInfo(sarks, '24', filters, thresholds)
#'
#' @export
blockInfo = function(sarks, block, filters, thresholds, kMax=12L) {
    if (length(block) > 1) {
        return(do.call(
            rbind,
            lapply(block, blockInfo, sarks=sarks,
                   filters=filters, thresholds=thresholds, kMax=kMax)
        ))
    }
    kMax = as.integer(kMax)
    javaOut = sarks$printBlockInfo(filters, thresholds$.java,
                                   block, NULL, kMax)
    out = utils::read.table(
        textConnection(javaOut),
        sep='\t', header=TRUE, row.names=NULL, na.strings='null',
        stringsAsFactors=FALSE, check.names=FALSE
    )
    return(out)
}


#' Slurp fasta file into character vector
#' 
#' Convenience function for quickly reading a fasta file into an R character vector.
#'
#' @param f the path for the fasta file to read in
#'
#' @return a named character vector: names are the sequence names from
#'     the fasta file, values are the sequences themselves.
#'
#' @export
slurpFasta = function(f) {
    ## tmpf = NULL
    ## if (grepl('https://raw\\.githubusercontent\\.com', f)) {
    ##     tmpf = tempfile()
    ##     utils::download.file(f, destfile=tmpf)
    ##     f = tmpf
    ## }
    slurped = paste(gsub('^>(.*)$', '>>>\\1>>>', readLines(f)), collapse='')
    ## if (length(tmpf) > 0) {file.remove(tmpf);}
    slurped = gsub('^>>>', '', slurped)
    lines = strsplit(slurped, '>>>')[[1]]
    return(structure(
        gsub('\\s+', '', lines[2*(1:(length(lines)/2))]),
        names = lines[-1 + 2*(1:(length(lines)/2))]
    ))
}


#' Reverse-complement DNA sequence
#' 
#' Takes character vector input representing DNA sequence and returns
#' reverse-complement as character vector.
#'
#' @param s character vector DNA sequence(s)
#'
#' @return character vector containing reverse-complement of s
#'
#' @export
dnaRevComp = function(s) {
    if (length(s) > 1) {return(sapply(s, dnaRevComp))}
    comp = c(a='t', A='T', c='g', C='G', g='c', G='C',
             n='N', N='N', t='a', T='A', `[`=']', `]`='[', `-`='-')
    return(paste(rev(comp[ strsplit(s, split='')[[1]] ]),
                 collapse = ''))
}


#' Count occurrences of regular expression
#'
#' Counts how often a regular expression (or vector of regular
#' expressions) occurs in each element of a character vector.
#'
#' @param regex character vector of regular expressions to search for
#'
#' @param seqs character vector of sequences in which to search for
#'     and count occurrences of regex
#'
#' @param overlap logical value: should overlapping occurrences of
#'     regex be counted as multiple hits?
#'
#' @return if length(regex) is one, returns integer vector of counts;
#'     if length(regex) is more than one, returns matrix of counts:
#'     one row per sequence in seqs, one column per expression in
#'     regex
#'
#' @export
regexCounts = function(regex, seqs, overlap=FALSE) {
    if (length(regex) > 1) {
        return(sapply(regex, regexCounts, seqs=seqs, overlap=overlap))
    }
    if (overlap) {regex = paste0('(?=', regex, ')')}
    return(nchar(gsub('[^_]+', '',
                      gsub(regex, '_', toupper(seqs), perl=TRUE))))
}


#' Count occurrences of k-mer
#'
#' Counts how often a k-mer (or vector of k-mers) occurs in each
#' element of a character vector.
#'
#' @param kmer character vector of k-mers to search for.
#'
#' @param seqs character vector of sequences in which to search for
#'     and count occurrences of kmer.
#'
#' @param directional logical value: if FALSE, counts occurrences of
#'     either kmer or its reverse-complement.
#'
#' @param overlap logical value: should overlapping occurrences of
#'     kmer be counted as multiple hits?
#'
#' @return if length(kmer) is one, returns integer vector of counts;
#'     if length(kmer) is more than one, returns matrix of counts:
#'     one row per sequence in seqs, one column per expression in
#'     regex
#'
#' @export
kmerCounts = function(kmer, seqs, directional=TRUE, overlap=FALSE) {
    if (is.numeric(kmer)) {
        nucs = c('A', 'C', 'G', 'T')
        nucs = lapply(1:kmer, function(...) {nucs})
        kmer = apply(X = do.call(expand.grid, nucs),
                     MARGIN = 1,
                     FUN = function(x) paste(x, collapse=''))
        names(kmer) = kmer
    }
    if (!directional) {
        kmer = structure(paste0(kmer, '|', dnaRevComp(kmer)), names=kmer)
    } else if (length(names(kmer)) == 0) {
        kmer = structure(kmer, names=kmer)
    }
    return(regexCounts(kmer, seqs, overlap))    
}


#' Count occurrences of k-mer clusters
#'
#' Counts how often any k-mer from a cluster of k-mers (or list of
#' clusters of k-mers) occurs in each element of a character vector.
#'
#' @param kmers character vector of k-mers composing cluster to search
#'     for, or a named list of such character vectors to count
#'     multiple clusters.
#'
#' @param seqs character vector of sequences in which to search for
#'     and count occurrences of kmers.
#'
#' @param directional logical value: if FALSE, counts occurrences of
#'     either cluster(s) of k-mers or their reverse-complements.
#'
#' @param overlap logical value: should overlapping occurrences of
#'     k-mers be counted as multiple hits?
#'
#' @return if cluster is a single character vector (of any length),
#'     returns integer vector of counts; if cluster is a list of
#'     character vectors, returns matrix of counts: one row per
#'     sequence in seqs, one column per character vector in cluster
#'
#' @export
clusterCounts = function(kmers, seqs, directional=TRUE, overlap=FALSE) {
    if (is.list(kmers)) {
        clCounts = lapply(kmers, clusterCounts,
                          seqs=seqs, directional=directional, overlap=overlap)
        out = matrix(0, nrow=length(seqs), ncol=length(kmers))
        rownames(out) = names(seqs)
        for (i in 1:length(kmers)) {
            out[names(clCounts[[i]]), i] = clCounts[[i]]
        }
        colnames(out) = names(kmers)
        return(out)
    }
    if (!directional) {
        kmers = sort(unique(c(kmers, sapply(kmers, dnaRevComp))))
    }
    pattern = paste(kmers, collapse='|')
    return(regexCounts(pattern, seqs, overlap))
}

#' Locate occurrences of regular expression
#'
#' Find locations of matches of a regular expression (or vector of
#' regular expressions) in each element of a character vector.
#'
#' @param regex character vector of regular expressions to search for
#'
#' @param seqs character vector of sequences in which to locate regex
#'
#' @return If only a single regex is searched for: data.frame with two
#'     columns: `seqid' containing the name of the sequence from seqs
#'     in which the regex was found and `location' giving the 1-based
#'     position at which the regex was found. If length(regex) greater
#'     than one, adds additional column `regex' indicating the name of
#'     the regex located.
#'
#' @export
regexLocate = function(regex, seqs) {
    if (length(regex) > 1) {
        if (length(names(regex)) == 0) {names(regex) = regex}
        out = lapply(regex, regexLocate, seqs=seqs)
        for (rn in names(regex)) {out[[rn]]$regex = rn}
        return(do.call(rbind, out)[ , c('seqid', 'regex', 'location')])
    }
    inseqs = seqs
    frags = strsplit(gsub(paste0('(', regex, ')'), '_\\1', toupper(seqs)), '_')
    for (i in 1:length(inseqs)) {
        names(frags[[i]]) = rep(names(inseqs)[i], length(frags[[i]]))
    }
    names(frags) = NULL
    fragLocs = sapply(sapply(frags, nchar), cumsum)
    fragLocs = unlist(sapply(fragLocs, function(.) {
        if (length(.) > 1) {
            return(.[1:(length(.)-1)] + 1)
        } else {
            return(integer(0))
        }
    }))
    return(data.frame(
        seqid = names(fragLocs),
        location = fragLocs
    ))
}


#' Locate occurrences of specified k-mers
#'
#' Find locations of matches of vector of k-mers in each element of a
#' character vector.
#'
#' @param kmers character vector of k-mers to search for
#'
#' @param seqs character vector of sequences in which to locate kmer
#'
#' @param directional logical value: if FALSE, counts occurrences of
#'     either kmers or their reverse-complements.
#'
#' @return data.frame with three columns: `seqid' containing the name
#'     of the sequence from seqs in which the k-mer was found; `kmer'
#'     indicating the k-mer located; and `location' giving the 1-based
#'     position at which the match was found.
#'
#' @export
locateKmers = function(kmers, seqs, directional=TRUE) {
    patterns = structure(kmers, names=kmers)
    if (!directional) {
        patterns = structure(paste0(kmers, '|', dnaRevComp(kmers)),
                             names=kmers)
    }
    locations = lapply(patterns, regexLocate, seqs=seqs)
    for (kmer in kmers) {locations[[kmer]]$kmer = kmer}
    out = do.call(rbind, locations)[ , c('seqid', 'kmer', 'location')]
    rownames(out) = NULL
    return(out)
}


#' Locate occurrences of specified clusters of k-mers
#'
#' Find locations of matches of list of character vectors of k-mers in
#' each element of a character vector.
#'
#' @param clusters list of character vectors of k-mers to search for
#'
#' @param seqs character vector of sequences in which to locate kmer
#'
#' @param directional logical value: if FALSE, counts occurrences of
#'     either k-mers within each cluster or their reverse-complements.
#'
#' @return data.frame with three columns: `seqid' containing the name
#'     of the sequence from seqs in which the match was found;
#'     `cluster' indicating the cluster from wich a k-mer was located;
#'     and `location' giving the 1-based position at which the match
#'     was found.
#'
#' @export
locateClusters = function(clusters, seqs, directional=TRUE) {
    out = list()
    for (cn in names(clusters)) {
        kmers = clusters[[cn]]
        if (!directional) {
            kmers = sort(unique(c(kmers, sapply(kmers, dnaRevComp))))
        }
        pattern = paste(kmers, collapse='|')
        locations = regexLocate(pattern, seqs)
        locations$cluster = cn
        out[[cn]] = locations
    }
    out = do.call(rbind, out)[ , c('seqid', 'cluster', 'location')]
    rownames(out) = NULL
    return(out)
}


## orientSeqs = function(seqs, representative, k=4) {
##     names(seqs) = seqs
##     seqKmers = kmerCounts(seqs, k, overlap=TRUE)
##     revComps = sapply(seqs, dnaRevComp)
##     revKmers = kmerCounts(revComps, k, overlap=TRUE)
##     reverse = rep(FALSE, length(seqs))
##     repKmers = as.numeric(seqKmers[representative, ])
##     for (i in 1:length(seqs)) {
##         if (sum(as.numeric(revKmers[i, ]) * repKmers) >
##             sum(as.numeric(seqKmers[i, ]) * repKmers)) {
##             reverse[i] = TRUE
##         }
##     }
##     return(ifelse(reverse, revComps, seqs))
## }


#' Cluster k-mers
#'
#' Takes a set of k-mer sequences and returns a list of partitioning
#' the input k-mers into clusters of more similar k-mers. Hierarchical
#' clustering (average linkage) is performed based on Jaccard
#' coefficient distance metric applied treating each k-mer as the set
#' of all tetramers which can be found as substrings within it.
#' Tetramers which are reverse-complements of each other are treated
#' as equivalent.
#'
#' @param kmers character vector of k-mers to partition into clusters
#'
#' @param k length of sub-k-mers (default k=4 to use tetramers) with
#'     which to calculate Jaccard distances for clustering
#'
#' @param nClusters number of clusters to partition kmers into; if set
#'     to NULL (default value), selects number of clusters to maximize
#'     the average silhouette score
#'     (\url{https://en.wikipedia.org/wiki/Silhouette_(clustering)}).
#'
#' @export
clusterKmers = function(kmers, k=4, nClusters=NULL) {
    kmers = unique(kmers)
    if (length(kmers) <= 3 &&
        (length(nClusters) == 0 || nClusters > 3)) {
        stop('Must specify nClusters <= 3 for length(kmers) <= 3.')
    }
    tetraCounts = kmerCounts(k, kmers, overlap=TRUE)
    unorientedTetramers = sort(unique(sapply(
        colnames(tetraCounts),
        function(km) {min(km, dnaRevComp(km))}
    )))
    tetraCountsUnorient = list()
    for (col in unorientedTetramers) {
        if (col == dnaRevComp(col)) {
            tetraCountsUnorient[[col]] = tetraCounts[ , col]
        } else {
            tetraCountsUnorient[[col]] = tetraCounts[ , col] +
                                         tetraCounts[ , dnaRevComp(col)]
        }
    }
    tetraCountsUnorient = data.frame(tetraCountsUnorient, row.names=kmers)
    tetraCountsUnorient = tetraCountsUnorient[ , colSums(tetraCountsUnorient) > 0]
    ## now calculate Jaccard distance matrix:
    pq = as.matrix(tetraCountsUnorient) %*%t(as.matrix(tetraCountsUnorient))
    p2 = rowSums(tetraCountsUnorient^2)
    d = 1 - pq / (outer(p2, p2, `+`) - pq)
    rownames(d) = colnames(d) = rownames(tetraCountsUnorient)
    diag(d) = Inf
    d[d == 0] = min(d[d > 0]) / 2
    diag(d) = 0
    ddist = stats::as.dist(d)
    ## cluster based on the distances in d
    hcout = stats::hclust(ddist, method='average')
    ## use silhouette method to determine number of clusters
    if (length(nClusters) == 0) {
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
        ctouts = data.frame(stats::cutree(hcout, k=2:(nrow(d)-1)), check.names=FALSE)
        nClusters = integralOptimize(
            f = function(nc) {
                mean(cluster:::silhouette.default(ctouts[ , as.character(nc)],
                                                  ddist)[ , 3])
            },
            lower = 2,
            upper = nrow(d)-1
        )
    }
    ctout = stats::cutree(hcout, k=nClusters)
    clusters = lapply(1:max(ctout), function(cl) {d[ctout==cl, ctout==cl, drop=FALSE]})
    clustReps = sapply(clusters, function(d) {rownames(d)[which.min(rowMeans(d))]})
    names(clusters) = clustReps
    return(lapply(clusters, function(d) {rownames(d)}))
}
