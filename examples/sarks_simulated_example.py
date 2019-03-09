#!/usr/bin/env python3

import numpy as np
import pandas as pd

from sarks import Sarks

## load scores into pandas Series
scores = pd.read_csv('simulated_scores.tsv', sep='\t', header=0, index_col=0).iloc[:, 0]
## initialize Sarks object
sarks = Sarks(
    inFasta = 'simulated_seqs.fa',          ## input sequences;
                                            ##   names should match scores.index
    catFasta = 'simulated_seqs_concat.fa',  ## name for intermediate
                                            ##   concatenated fasta file
    scores = scores,
    halfWindow = 4                          ## kernel half-width kappa
)

topLocations = sarks.peaks(theta = 1,
                           minGini = 0)
topTable = sarks.subtable(topLocations).sort_values('khat',
                                                    ascending=False)
#    i     s        kmer    khat block   wi      gini  score  windowed
                                                                    
# 2257  3959  CATACTGAGA  10.250    22  194  0.888889      1       1.0
# 2258  4518  CATACTGAGA  10.250    25    0  0.888889      1       1.0
# 2256  3544  CATACTGAGA   9.625    21   30  0.864198      1       1.0
# 1460  3960   ATACTGAGA   9.250    22  195  0.888889      1       1.0
# 1461  4519   ATACTGAGA   9.250    25    1  0.888889      1       1.0
# 1459  3545   ATACTGAGA   8.750    21   31  0.888889      1       1.0
# 1462  3456    ATACTGAG   8.500    20  193  0.864198      1       1.0
# 1458  4442    ATACTGAG   8.250    24  175  0.864198      1       1.0
# 5864  3961    TACTGAGA   8.250    22  196  0.888889      1       1.0
# 5865  4520    TACTGAGA   8.250    25    2  0.888889      1       1.0
# 1463  5595    ATACTGAG   7.875    29   73  0.864198      1       1.0
# 5863  3546    TACTGAGA   7.750    21   32  0.888889      1       1.0
# 5862  4443     TACTGAG   7.250    24  176  0.864198      1       1.0
# 1464  5174     ATACTGA   7.125    27  154  0.839506      1       1.0
# 5861  5430     TACTGAG   6.875    28  159  0.839506      1       1.0
# 1465  4232      ATACTG   6.250    23  216  0.814815      1       1.0

extTopTable = sarks.extendKmers(topTable)
#    i     s        kmer    khat block   wi      gini  score  windowed
                                                                    
# 2257  3959  CATACTGAGA  10.250    22  194  0.888889      1       1.0
# 2258  4518  CATACTGAGA  10.250    25    0  0.888889      1       1.0
# 2256  3544  CATACTGAGA   9.625    21   30  0.864198      1       1.0
# 2257  3959  CATACTGAGA     NaN    22  194       NaN      1       NaN
# 2258  4518  CATACTGAGA     NaN    25    0       NaN      1       NaN
# 2256  3544  CATACTGAGA     NaN    21   30       NaN      1       NaN
# 2259  3455  CATACTGAGA     NaN    20  192       NaN      1       NaN
# 2255  4441  CATACTGAGA     NaN    24  174       NaN      1       NaN
# 2257  3959  CATACTGAGA     NaN    22  194       NaN      1       NaN
# 2258  4518  CATACTGAGA     NaN    25    0       NaN      1       NaN
# 2260  5594  CATACTGAGA     NaN    29   72       NaN      1       NaN
# 2256  3544  CATACTGAGA     NaN    21   30       NaN      1       NaN
# 2255  4441  CATACTGAGA     NaN    24  174       NaN      1       NaN
# 2261  5173  CATACTGAGA     NaN    27  153       NaN      1       NaN
# 2254  5428  CATACTGAGA     NaN    28  157       NaN      1       NaN
# 2262  4231  CATACTGAGA     NaN    23  215       NaN      1       NaN
