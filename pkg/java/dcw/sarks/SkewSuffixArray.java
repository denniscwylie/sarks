package dcw.sarks;

/**
 * Implementation of skew algorithm for construction of suffix array
 * based on code provided in: "Kärkkäinen, J., & Sanders, P. (2003,
 * June). Simple linear work suffix array construction. In
 * International Colloquium on Automata, Languages, and Programming
 * (pp. 943-955). Springer, Berlin, Heidelberg."
 */
public class SkewSuffixArray {
    private static boolean leq(int a1, int a2, int b1, int b2) {
        return (a1 < b1 || (a1 == b1 && a2 <= b2));
    }
    private static boolean leq(int a1, int a2, int a3, int b1, int b2, int b3) {
        return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
    }

    /**
     * stably sort indexes from a[0..n-1] to b[0..n-1] with keys in 0..K from r
     * (indices of r shifted by rShift positions)
     */
    private static void radixPass(int[] a, int[] b, int[] r, int rShift, int n, int K) {
        // count occurrences
        int[] c = new int[K + 1];
        // reset counters
        for (int i = 0; i <= K; i++) {c[i] = 0;}
        // count occurrences
        for (int i = 0; i < n; i++) {c[r[rShift + a[i]]]++;}
        // exclusive prefix sums
        for (int i = 0, sum = 0; i <= K; i++) {
            int t = c[i]; c[i] = sum; sum += t;
        }
        // sort
        for (int i = 0; i < n; i++) {b[c[r[rShift + a[i]]]++] = a[i];}
    }

    /**
     * find the suffix array SA of s[0..n-1] in {1..K}^n
     * require s[n]=s[n+1]=s[n+2]=0, n >= 2.
     */
    static void suffixArray(int[] s, int[] SA, int n, int K) {
        int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
        int[] s12 = new int[n02 + 3]; s12[n02] = s12[n02+1] = s12[n02+2] = 0;
        int[] SA12 = new int[n02 + 3]; SA12[n02] = SA12[n02+1] = SA12[n02+2] = 0;
        int[] s0 = new int[n0];
        int[] SA0 = new int[n0];
        
        // generate positions of mod 1 and mod 2 suffixes
        // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
        for (int i=0, j=0; i<n+(n0-n1); i++) {
            if ((i%3) != 0) {s12[j++] = i;}
        }

        // lsb radix sort the mod 1 and mod 2 triples
        radixPass(s12, SA12, s, 2, n02, K);
        radixPass(SA12, s12, s, 1, n02, K);
        radixPass(s12, SA12, s, 0, n02, K);

        // find lexicographic names of triples
        int name = 0, c0 = -1, c1 = -1, c2 = -1;
        for (int i = 0; i < n02; i++) {
            if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
                name++; c0 = s[SA12[i]]; c1 = s[SA12[i]+1]; c2 = s[SA12[i]+2];
            }
            if ((SA12[i] % 3) == 1) {
                // left half
                s12[SA12[i] / 3] = name;
            }
            else {
                // right half
                s12[SA12[i] / 3 + n0] = name;
            }
        }

        // recurse if names are not yet unique
        if (name < n02) {
            suffixArray(s12, SA12, n02, name);
            // store unique names in s12 using the suffix array
            for (int i = 0; i < n02; i++) {s12[SA12[i]] = i + 1;}
        } else {
            // generate the suffix array of s12 directly
            for (int i = 0; i < n02; i++) {SA12[s12[i] - 1] = i;}
        }

        // stably sort the mod 0 suffixes from SA12 by their first character
        for (int i=0, j=0; i<n02; i++) {
            if (SA12[i] < n0) {s0[j++] = 3 * SA12[i];}
        }
        radixPass(s0, SA0, s, 0, n0, K);

        // merge sorted SA0 suffixes and sorted SA12 suffixes
        for (int p=0, t=n0-n1, k=0; k<n; k++) {
            // pos of current offset 12 suffix
            int i = (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2);
            // pos of current offset 0 suffix
            int j = SA0[p];

            if (SA12[t] < n0 ?
                leq(s[i], s12[SA12[t]+n0], s[j], s12[j/3]) :
                leq(s[i], s[i+1], s12[SA12[t]-n0+1], s[j], s[j+1], s12[j/3+n0])) {
                // suffix from SA12 is smaller
                SA[k] = i; t++;
                if (t == n02) {
                    // done --- only SA0 suffixes left
                    for (k++; p < n0; p++, k++) {SA[k] = SA0[p];}
                }
            } else {
                SA[k] = j; p++;
                if (p == n0) {
                    // done --- only SA12 suffixes left
                    for (k++; t < n02; t++, k++) {
                        SA[k] = (SA12[t] < n0 ?
                                 SA12[t] * 3 + 1 :
                                 (SA12[t] - n0) * 3 + 2);
                    }
                }
            }
        }
    }
}
