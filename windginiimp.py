#!/usr/bin/env python3

import ctypes
import inline
import numpy as np
import pandas as pd
from numpy.ctypeslib import ndpointer
from scipy import stats


## -----------------------------------------------------------------
windGiniImpC = inline.c("""
#include <math.h>

int* blockCounts(int i, int halfWindow, int* block, int nBlocks) {
    int* out = (int*)calloc(nBlocks, sizeof(int));
    int j;
    for (j=0; j<nBlocks; j++) {out[j] = 0;}
    for (j=-halfWindow; j<=halfWindow; j++) {
        out[block[i+j]]++;
    }
    return out;
}

double giniImpurity(int* counts, int nBlocks) {
    double out = 0;
    int total = 0;
    int b;
    for (b=0; b<nBlocks; b++) {
        total += counts[b];
    }
    for (b=0; b<nBlocks; b++) {
        out += (double)counts[b] * (total - counts[b]);
    }
    out /= (double)(total * total);
    return out;
}

double* giniImpurities(int halfWindow, int* block, int lenBlocks, int nBlocks) {
    int* bCounts = blockCounts(halfWindow, halfWindow, block, nBlocks);
    int halfWindow2 = 2 * halfWindow;
    int nGinis = lenBlocks - halfWindow2;
    double* ginis = (double*)calloc(nGinis, sizeof(double));
    ginis[0] = giniImpurity(bCounts, nBlocks);
    int total = (2 * halfWindow) + 1;
    double total2 = (double)(total * total);
    int dGini = 0;
    int old, new, i;
    for (i=1; i<nGinis; i++) {
        old = block[i-1];
        new = block[i+halfWindow2];
        dGini = -bCounts[old] * (total - bCounts[old]);
        if (bCounts[new] > 0) {
            dGini -= bCounts[new] * (total - bCounts[new]);
        }
        bCounts[old]--;
        bCounts[new]++;
        if (bCounts[old] > 0) {
            dGini += bCounts[old] * (total - bCounts[old]);
        }
        dGini += bCounts[new] * (total - bCounts[new]);
        ginis[i] = ginis[i-1] + ((double)dGini / total2);
    }
    return ginis;
}
""")


def windGiniImpurities(b, halfWindow):
    bInt = pd.Series(range(len(b.unique())),
                     index = b.unique()).loc[b]
    nGinis = len(bInt) - 2*halfWindow
    windGiniImpC.giniImpurities.restype = ndpointer(dtype = ctypes.c_double,
                                                    shape = (nGinis,))
    bType = ctypes.c_int * len(bInt)
    ginis = windGiniImpC.giniImpurities(halfWindow,
                                        bType(*bInt),
                                        len(bInt),
                                        len(bInt.unique()))
    ginis = pd.Series(np.array(ginis),
                      index = b.index[halfWindow:(len(b)-halfWindow)])
    return ginis
