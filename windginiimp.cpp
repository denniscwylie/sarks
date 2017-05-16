#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using std::ios;
using std::cout;

int* blockCounts(int i, int halfWindow, std::vector<int> block, int nBlocks) {
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

double* giniImpurities(int halfWindow, std::vector<int> block) {
    int lenBlocks = block.size();
    std::set<int> blockSet(block.begin(), block.end());
    int nBlocks = blockSet.size();
    int* bCounts = blockCounts(halfWindow, halfWindow, block, nBlocks);
    int halfWindow2 = 2 * halfWindow;
    int nGinis = lenBlocks - halfWindow2;
    double* ginis = (double*)calloc(nGinis, sizeof(double));
    ginis[0] = giniImpurity(bCounts, nBlocks);
    int total = (2 * halfWindow) + 1;
    double total2 = (double)(total * total);
    int dGini = 0;
    int oldb, newb, i;
    for (i=1; i<nGinis; i++) {
        oldb = block[i-1];
        newb = block[i+halfWindow2];
        dGini = -bCounts[oldb] * (total - bCounts[oldb]);
        if (bCounts[newb] > 0) {
            dGini -= bCounts[newb] * (total - bCounts[newb]);
        }
        bCounts[oldb]--;
        bCounts[newb]++;
        if (bCounts[oldb] > 0) {
            dGini += bCounts[oldb] * (total - bCounts[oldb]);
        }
        dGini += bCounts[newb] * (total - bCounts[newb]);
        ginis[i] = ginis[i-1] + ((double)dGini / total2);
    }
    return ginis;
}

int main(int argc, char** argv) {
    std::fstream in(argv[2], ios::in);
    std::vector<int> blocks;
    for (int block; in >> block;) {
        blocks.push_back(block);
    }
    int halfWindow = std::atoi(argv[1]);
    int nGinis = blocks.size() - (2 * halfWindow);
    double* giniImp = giniImpurities(halfWindow, blocks);
    for (int i=0; i<nGinis; i++) {
        cout << std::setprecision(10) << giniImp[i] << "\n";
    }
    return 0;
}
