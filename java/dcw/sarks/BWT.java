package dcw.sarks;

import java.util.HashMap;
import java.util.TreeMap;

/**
 * Implementation of Burrows-Wheeler transform for purpose of exact
 * k-mer matching based on backward search algorithm described in:
 * "Ferragina, P., & Manzini, G. (2000, November). Opportunistic data
 * structures with applications. In Proceedings 41st Annual Symposium
 * on Foundations of Computer Science (pp. 390-398). IEEE."
 */
public class BWT {

    private String seq;
    private int[] sa;
    private char[] bwt;
    private TreeMap<Character,Integer> c;
    private HashMap<Character,int[]> o;

    private static final int BLOCK_SIZE = 100;

    public BWT(String seq0, int[] sa0) {
        this.seq = seq0;
        this.sa = sa0;
        this.bwt = new char[this.seq.length()];
        this.c = new TreeMap<Character,Integer>();
        this.o = new HashMap<Character,int[]>();
        for (int i=0; i<this.bwt.length; i++) {
            if (this.sa[i] > 0) {
                this.bwt[i] = this.seq.charAt(this.sa[i] - 1);
            } else {
                this.bwt[i] = '$';
            }
            char ichar = this.bwt[i];
            if (!this.c.containsKey(ichar)) {this.c.put(ichar, 0);}
            this.c.put(ichar, this.c.get(ichar)+1);
            if (!this.o.containsKey(ichar)) {
                int[] oPart = new int[((this.bwt.length-1) / BLOCK_SIZE)+1];
                for (int j=0; j<oPart.length; j++) {
                    oPart[j] = 0;
                }
                this.o.put(ichar, oPart);
            }
            this.o.get(ichar)[i / BLOCK_SIZE]++;
        }
        int cumul = 0;
        int oSize = ((this.bwt.length-1) / BLOCK_SIZE) + 1;
        for (Character ichar : this.c.keySet()) {
            cumul += this.c.get(ichar);
            this.c.put(ichar, cumul - this.c.get(ichar));
            int[] ochar = this.o.get(ichar);
            for (int i=1; i<oSize; i++) {ochar[i] += ochar[i-1];}
        }
    }

    public int[] findKmer(String kmer) {
        int k = kmer.length();
        int start = 0;
        int end = this.seq.length();
        for (int i=k-1; i>=0; i--) {
            char ichar = kmer.charAt(i);
            int startRounded = (start / BLOCK_SIZE) - 1;
            int oStartRounded = (startRounded < 0 ? 0 : this.o.get(ichar)[startRounded]);
            int newStart = this.c.get(ichar) + oStartRounded;
            for (int j=(BLOCK_SIZE*(startRounded+1)); j<start; j++) {
                if (this.bwt[j] == ichar) {newStart++;}
            }
            start = newStart;
            int endRounded = (end / BLOCK_SIZE) - 1;
            int oEndRounded = (endRounded < 0 ? 0 : this.o.get(ichar)[endRounded]);
            int newEnd = this.c.get(ichar) + oEndRounded;
            for (int j=(BLOCK_SIZE*(endRounded+1)); j<end; j++) {
                if (this.bwt[j] == ichar) {newEnd++;}
            }
            end = newEnd;
        }
        return new int[] {start, end};
    }
}
