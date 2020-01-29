package dcw.sarks;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class Sarks {

    private int nThreads;

    private int halfWindow;
    private int spatialLength;
        
    private double[] scores;
    private String[] transcripts;
    private HashMap<String,Integer> transcriptPosition;
    private String catSeq;
    private int[] bounds;
    private int[] sa;
    private int[] saInv;
    private float[] windGini;
    private float[] spatGini;
    private float[] windowed;
    private float[] spatialWindowed;

    private BWT bwt;

    public Sarks(String inFasta, String scoreFile,
                 int halfWindow0, int spatialLength0,
                 int nThreads0, boolean burrowsWheeler) {
        this.nThreads = nThreads0;
        this.halfWindow = halfWindow0;
        this.spatialLength = spatialLength0;
        HashMap<String,String> seqs = null;
        TreeMap<String,Double> scoreMap = null;
        try {
            seqs = SarksUtilities.readFasta(inFasta);
            scoreMap = SarksUtilities.readScores(scoreFile);
        } catch (Exception e) {}
        scoreMap.keySet().retainAll(seqs.keySet());
        seqs.keySet().retainAll(scoreMap.keySet());
        this.transcripts = new String[scoreMap.size()];
        this.transcriptPosition = new HashMap<String,Integer>();
        this.scores = new double[scoreMap.size()];
        int i = 0;
        for (Map.Entry<String,Double> entry : scoreMap.entrySet()) {
            this.transcripts[i] = entry.getKey();
            this.transcriptPosition.put(this.transcripts[i], i);
            this.scores[i] = entry.getValue();
            i++;
        }
        this.concatenateSeqs(seqs);
        this.calcSuffixArray();
        if (burrowsWheeler) {
            this.bwt = new BWT(this.catSeq, this.sa);
        } else {
            this.bwt = null;
        }
        this.windowGini();
        this.window();
        if (this.spatialLength > 1) {
            this.spatialGini();
            this.spatialWindow();
        } else {
            this.spatGini = null;
            this.spatialWindowed = null;
        }
    }
    public Sarks(String inFasta, String scoreFile, int halfWindow0) {
        this(inFasta, scoreFile, halfWindow0, 0, 1, false);
    }

    public Sarks(Sarks sarks, Integer[] permutation) {
        this.nThreads = sarks.nThreads;
        this.halfWindow = sarks.halfWindow;
        this.spatialLength = sarks.spatialLength;
        this.transcripts = sarks.transcripts;
        this.transcriptPosition = sarks.transcriptPosition;
        this.scores = new double[sarks.scores.length];
        for (int b=0; b<scores.length; b++) {
            this.scores[ permutation[b] ] = sarks.scores[b];
        }
        this.catSeq = sarks.catSeq;
        this.bounds = sarks.bounds;
        this.sa = sarks.sa;
        this.saInv = sarks.saInv;
        this.bwt = sarks.bwt;
        this.windGini = sarks.windGini;
        this.spatGini = sarks.spatGini;
        this.window();
        if (this.spatialLength > 1) {
            this.spatialWindow();
        }
    }

    private void concatenateSeqs(HashMap<String,String> seqs) {
        StringBuilder sb = new StringBuilder();
        this.bounds = new int[this.transcripts.length + 1];
        this.bounds[0] = 0;
        int i = 1;
        for (String t : this.transcripts) {
            String tSeq = seqs.get(t);
            sb.append(tSeq);
            sb.append("$");
            this.bounds[i] = this.bounds[i-1] + tSeq.length() + 1;
            i++;
        }
        this.catSeq = sb.toString();        
    }

    private void calcSuffixArray() {
        int[] chars = new int [this.catSeq.length() + 3];
        for (int i=this.catSeq.length()-1; i>=0; i--) {
            chars[i] = this.catSeq.charAt(i);
        }
        for (int i=this.catSeq.length(); i<chars.length; i++) {chars[i] = 0;}
        int alphabetSize = SarksUtilities.max(chars, this.catSeq.length());
        int[] SA = new int [chars.length];
        dcw.sarks.SkewSuffixArray.suffixArray(
                chars, SA, this.catSeq.length(), alphabetSize);
        this.sa = new int[this.catSeq.length()];
        this.saInv = new int[this.sa.length];
        for (int i=0; i<this.sa.length; i++) {
            this.sa[i] = SA[i];
            this.saInv[ SA[i] ] = i;
        }
    }

    public void reset(int halfWindow0, Integer spatialLength0, boolean smooth) {
        this.halfWindow = halfWindow0;
        if (spatialLength0 != null) {
            this.spatialLength = spatialLength0;
        }
        this.windowGini();
        if (smooth) {
            this.window();
        } else {
            this.windowed = null;
        }
        if (this.spatialLength > 1) {
            this.spatialGini();
            if (smooth) {
                this.spatialWindow();
            } else {
                this.spatialWindowed = null;
            }
        } else {
            this.spatGini = null;
            this.spatialWindowed = null;
        }        
    }
    public void reset(int halfWindow0) {this.reset(halfWindow0, null, true);}

    public void resetSpatial(int spatialLength0, boolean smooth) {
        this.spatialLength = spatialLength0;
        this.spatialGini();
        if (smooth) {
            this.spatialWindow();
        } else {
            this.spatialWindowed = null;
        }
    }
    public void resetSpatial(int spatialLength0) {this.resetSpatial(spatialLength0, true);}


    // // -------------------------------------------------------------------------
    // public int[] invertedSuffixArray() {
    //     int[] saInv = new int[this.sa.length];
    //     for (int i=0; i<this.sa.length; i++) {
    //         saInv[ this.sa[i] ] = i;
    //     }
    //     return saInv;
    // }


    // -------------------------------------------------------------------------
    public int sourceBlock(int s, int mindex, int maxdex) {
        if (maxdex == (mindex + 1)) {return mindex;}
        int middex = (mindex + maxdex) / 2;
        // middex++;
        if (s < this.bounds[middex]) {
            return this.sourceBlock(s, mindex, middex);
        } else {
            return this.sourceBlock(s, middex, maxdex);
        }
    }
    public int sourceBlock(int s) {return(this.sourceBlock(s, 0, this.bounds.length));}
    public int[] sourceBlock(int[] s0) {
        if (s0 != null) {
            int[] b = new int[s0.length];
            for (int i=0; i<s0.length; i++) {b[i] = this.sourceBlock(s0[i]);}
            return b;
        } else {
            int[] b = new int[this.sa.length];
            for (int block=0; block<this.transcripts.length; block++) {
                int nextBlock = block + 1;
                for (int s=this.bounds[block]; s<this.bounds[nextBlock]; s++) {
                    b[ this.saInv[s] ] = block;
                }
            }
            return b;
        }
    }

    public double[] sourceScore(int[] s) {
        int[] srcBlk = this.sourceBlock(s);
        double[] out = new double[srcBlk.length];
        for (int i=0; i<srcBlk.length; i++) {
            out[i] = this.scores[ srcBlk[i] ];
        }
        return out;
    }


    // -------------------------------------------------------------------------
    public int[] blockCounts(int i, int[] block) {
        int nBlocks = this.transcripts.length;
        int[] out = new int[nBlocks];
        for (int j=0; j<nBlocks; j++) {out[j] = 0;}
        for (int j=-this.halfWindow; j<=this.halfWindow; j++) {
            out[ block[i+j] ]++;
        }
        return out;
    }
    
    public static double giniImpurity(int[] counts) {
        double out = 0;
        int total = 0;
        for (int b=0; b<counts.length; b++) {
            total += counts[b];
        }
        for (int b=0; b<counts.length; b++) {
            out += (double)counts[b] * (total - counts[b]);
        }
        out /= (double)(total * total);
        return out;
    }

    private void windowGini() {
        int[] block = this.sourceBlock(null);
        this.windGini = new float[this.sa.length];
        int[] bCounts = this.blockCounts(this.halfWindow, block);
        double runningGini = this.giniImpurity(bCounts);
        this.windGini[this.halfWindow] = (float)runningGini;
        int total = (2 * this.halfWindow) + 1;
        double total2 = (double)(total * total);
        int dGini = 0;
        int oldb, newb;
        for (int i=(this.halfWindow+1); i<(this.sa.length-this.halfWindow); i++) {
            oldb = block[i-1-this.halfWindow];
            newb = block[i+this.halfWindow];
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
            runningGini += ((double)dGini / total2);
            this.windGini[i] = (float)runningGini;
        }
        for (int i=0; i<this.halfWindow; i++) {
            this.windGini[i] = this.windGini[this.halfWindow];
            this.windGini[this.sa.length-1-i] =
                    this.windGini[this.sa.length-1-this.halfWindow];
        }
    }

    private void spatialGini() {
        double windSize = (double)this.spatialLength;
        this.spatGini = new float[this.windGini.length];
        double runningMean = 0;
        for (int s=0; s<this.spatialLength; s++) {
            runningMean += this.windGini[ this.saInv[s] ];
        }
        runningMean /= windSize;
        this.spatGini[0] = (float)runningMean;
        for (int s=1; s<(this.sa.length-this.spatialLength); s++) {
            runningMean -= (this.windGini[ this.saInv[s-1] ] / windSize);
            runningMean += (this.windGini[ this.saInv[s+this.spatialLength] ] / windSize);
            this.spatGini[ this.saInv[s] ] = (float)runningMean;
        }
        for (int s=(this.saInv.length-this.spatialLength); s<this.sa.length; s++) {
            this.spatGini[s] = (float)runningMean;
        }
    }


    // -------------------------------------------------------------------------
    private void window() {
        double[] saScores = this.sourceScore(null);
        double windSize = (double)((2 * this.halfWindow) + 1);
        this.windowed = new float[this.sa.length];
        double runningMean = 0;
        for (int i=0; i<((2*this.halfWindow)+1); i++) {
            runningMean += saScores[i];
        }
        runningMean /= windSize;
        this.windowed[this.halfWindow] = (float)runningMean;
        for (int i=(this.halfWindow+1); i<(this.sa.length-this.halfWindow); i++) {
            runningMean -= (saScores[i-1-this.halfWindow] / windSize);
            runningMean += (saScores[i+this.halfWindow] / windSize);
            this.windowed[i] = (float)runningMean;
        }
        for (int i=0; i<this.halfWindow; i++) {
            this.windowed[i] = this.windowed[this.halfWindow];
            this.windowed[this.sa.length-1-i] = (float)runningMean;
        }
    }

    private void spatialWindow() {
        double windSize = (double)this.spatialLength;
        this.spatialWindowed = new float[this.sa.length];
        double runningMean = 0;
        for (int s=0; s<this.spatialLength; s++) {
            runningMean += this.windowed[ this.saInv[s] ];
        }
        runningMean /= windSize;
        this.spatialWindowed[0] = (float)runningMean;
        for (int s=1; s<(this.sa.length-this.spatialLength); s++) {
            runningMean -= (this.windowed[ this.saInv[s-1] ] / windSize);
            runningMean += (this.windowed[ this.saInv[s+this.spatialLength-1] ] / windSize);
            this.spatialWindowed[ this.saInv[s] ] = (float)runningMean;
        }
        for (int s=(this.sa.length-this.spatialLength); s<this.sa.length; s++) {
            this.spatialWindowed[ this.saInv[s] ] = (float)runningMean;
        }
    }


    // -------------------------------------------------------------------------
    private boolean startGood(int start, String kmer) {
        int k = kmer.length();
        int s = this.sa[start];
        int spk = s + k;
        if (spk > this.catSeq.length()) {spk = this.catSeq.length();}
        String kms = this.catSeq.substring(s, spk);
        if (kms.equals(kmer) && (start == 0)) {
            return true;
        } else if (kms.compareTo(kmer) >= 0) {
            int sm1 = this.sa[start-1];
            int sm1pk = sm1 + k;
            if (sm1pk > this.catSeq.length()) {sm1pk = this.catSeq.length();}
            String kmsm1 = this.catSeq.substring(sm1, sm1pk);
            if (kmsm1.compareTo(kmer) < 0) {
                return true;
            }
        }
        return false;
    }

    private boolean endGood(int end, String kmer) {
        int k = kmer.length();
        int sm1 = this.sa[end-1];
        int sm1pk = sm1 + k;
        if (sm1pk > this.catSeq.length()) {sm1pk = this.catSeq.length();}
        String kmem1 = this.catSeq.substring(sm1, sm1pk);
        if (kmem1.equals(kmer) && (end == this.sa.length)) {
            return true;
        } else if (kmem1.compareTo(kmer) <= 0) {
            int s = this.sa[end];
            int spk = s + k;
            if (spk > this.catSeq.length()) {spk = this.catSeq.length();}
            String kme = this.catSeq.substring(s, spk);
            if (kme.compareTo(kmer) > 0) {
                return true;
            }
        }
        return false;
    }
    
    public int[] findKmer(String kmer, int[] limits) {
        if (kmer.equals("")) {return new int[] {0, this.sa.length};}
        if (this.bwt != null) {return this.bwt.findKmer(kmer);}
        int k = kmer.length();
        int low = 0;
        int high = this.sa.length;
        if (limits != null) {
            low = limits[0];
            high = limits[1];
        }
        int start = low + (high-low) / 2;
        String km = null;
        int sstartpk;
        while (start < high && !this.startGood(start, kmer)) {
            sstartpk = this.sa[start]+k;
            if (sstartpk > this.catSeq.length()) {
                sstartpk = this.catSeq.length();
            }
            km = this.catSeq.substring(this.sa[start], sstartpk);
            if (km.compareTo(kmer) >= 0) {
                high = start;
            } else {
                low = start;
            }
            int newStart = low + (high - low) / 2;
            if (newStart == start) {
                start = high;
            } else {
                start = newStart;
            }
        }
        low = start;
        high = this.sa.length;
        int end = low + (high - low) / 2;
        int sendpk;
        while (end > start && !this.endGood(end, kmer)) {
            sendpk = this.sa[end] + k;
            if (sendpk > this.catSeq.length()) {
                sendpk = this.catSeq.length();
            }
            km = this.catSeq.substring(this.sa[end], sendpk);
            if (end == this.sa.length) {
                high = end;
            } else if (km.compareTo(kmer) <= 0) {
                low = end;
            } else {
                high = end;
            }
            int newEnd = low + (high - low) / 2;
            if (newEnd == end) {
                end = high;
            } else {
                end = newEnd;
            }
        }
        return new int[] {start, end};
    }
    public int[] findKmer(String kmer) {return this.findKmer(kmer, null);}

    public String kmer(int s, int k, int k0) {
        int catSeqMaxPos = this.catSeq.length() - 1;
        int s1 = s + k0; int s2 = s + k;
        if (s1 < 0) {s1 = 0;}
        if (s1 > catSeqMaxPos-1) {s1 = catSeqMaxPos-1;}
        if (s2 < 0) {s2 = 0;}
        if (s2 > catSeqMaxPos) {s2 = catSeqMaxPos;}
        return this.catSeq.substring(s1, s2);
    }
    public String kmer(int s, int k) {return this.kmer(s, k, 0);}
    
    public String[] kmers(int[] s, int k, int k0) {
        String[] out = new String[s.length];
        int catSeqMaxPos = this.catSeq.length() - 1;
        for (int sindex=0; sindex<s.length; sindex++) {
            int s1 = s[sindex] + k0;
            if (s1 < 0) {s1 = 0;}
            if (s1 > catSeqMaxPos-1) {s1 = catSeqMaxPos-1;}
            int s2 = s[sindex] + k;
            if (s2 < 0) {s2 = 0;}
            if (s2 > catSeqMaxPos) {s2 = catSeqMaxPos;}
            out[sindex] = this.catSeq.substring(s1, s2);
        }
        return out;
    }
    public String[] kmers(int[] s, int k) {return this.kmers(s, k, 0);}

    public double prefixAgreeSum(int i, int kmax) {
        String kmer = this.kmers(new int[] {this.sa[i]}, kmax)[0];
        int agreeSum = 0;
        int wstart = i - this.halfWindow;
        int wend = i + this.halfWindow + 1;
        int kstart = i;
        int kend = i+1;
        for (int k=kmax; k>0; k--) {
            int[] kwin = this.findKmer(kmer.substring(0, k));
            if (kwin[0] < wstart) {kwin[0] = wstart;}
            if (kwin[1] > wend) {kwin[1] = wend;}
            agreeSum += (k * ((kstart-kwin[0]) + (kwin[1]-kend)));
            kstart = kwin[0]; 
            kend = kwin[1];
        }
        return (double)agreeSum / (2.0 * this.halfWindow);
    }


    // -------------------------------------------------------------------------
    public ArrayList<Integer> filter(
            Float theta, Double minGini,
            Float spatialTheta, Double minSpatialGini,
            Boolean peakify) {
        ArrayList<Integer> pos = new ArrayList<Integer>();
        if (minGini != null) {
            if (minGini >= 1) {
                minGini = 1 - (1 - SarksUtilities.median(this.windGini)) * minGini;
            }
        } else {
            minGini = 0.;
        }
        if (theta == null) {theta = Float.NEGATIVE_INFINITY;}
        if (minSpatialGini != null) {
            if (minSpatialGini >= 1) {
                minSpatialGini =
                        1 - (1 - SarksUtilities.median(this.windGini)) * minSpatialGini;
            }
        } else {
            minSpatialGini = 0.;
        }
        if (peakify == null) {peakify = false;}
        boolean keep;
        int iLeft, iRight, s;
        for (int i=0; i<this.sa.length; i++) {
            if ((this.windowed[i] >= theta) && (this.windGini[i] >= minGini)) {
                keep = true;
                s = this.sa[i];
                if (spatialTheta != null) {
                    if ((this.spatialWindowed[i] < spatialTheta) ||
                        (this.spatGini[i] < minSpatialGini)) {
                        keep = false;
                    }
                }
                if (keep && peakify) {
                    iLeft = this.saInv[0];
                    if (s > 0) {
                        iLeft = this.saInv[s-1];
                    }
                    iRight = this.saInv[this.saInv.length-1];
                    if (s < (this.saInv.length-1)) {
                        iRight = this.saInv[s+1];
                    }
                    if ((this.windowed[iLeft] > this.windowed[i]) ||
                        (this.windowed[iRight] > this.windowed[i])) {
                        keep = false;
                    }
                }
                if (keep) {
                    pos.add(i);
                }
            }
        }
        return pos;
    }
    public ArrayList<Integer> filter(Float theta, Double minGini, Boolean peakify) {
        return this.filter(theta, minGini, null, null, peakify);
    }
    public ArrayList<Integer> filter(Float theta, Double minGini) {
        return this.filter(theta, minGini, null);
    }

    public ArrayList<ArrayList<Integer>> filter(ArrayList<HashMap> filters,
                                                Float[][] thresholds,
                                                Boolean peakify) {
        int hw0 = this.halfWindow;
        int sl0 = this.spatialLength;
        ArrayList<ArrayList<Integer>> out = new ArrayList<ArrayList<Integer>>();
        for (int f=0; f<filters.size(); f++) {
            HashMap filt = filters.get(f);
            Integer hw = (Integer)filt.get("halfWindow");
            Double ming = (Double)filt.get("minGini");
            Integer sl = (Integer)filt.get("spatialLength");
            Double minsg = (Double)filt.get("minSpatialGini");
            Float theta = null;
            Float spatialTheta = null;
            if (thresholds != null) {
                theta = thresholds[f][0];
                spatialTheta = thresholds[f][1];
            } else {
                theta = (Float)filt.get("theta");
                spatialTheta = (Float)filt.get("spatialTheta");
            }
            if ((hw != null) && (hw != this.halfWindow)) {
                this.reset(hw, sl, true);
            } else if ((sl != null) && (sl != this.spatialLength)) {
                this.resetSpatial(sl, true);
            }
            out.add(this.filter(theta, ming, spatialTheta, minsg, peakify));
        }
        if (this.halfWindow != hw0) {
            this.reset(hw0, sl0, true);
        } else if (this.spatialLength != sl0) {
            this.resetSpatial(sl0, true);
        }
        return out;
    }
    public ArrayList<ArrayList<Integer>> filter(ArrayList<HashMap> filters,
                                                Float[][] thresholds) {
        return this.filter(filters, thresholds, null);
    }
    public ArrayList<ArrayList<Integer>> filter(ArrayList<HashMap> filters,
                                                 Boolean peakify) {
        return this.filter(filters, null, peakify);
    }
    public ArrayList<ArrayList<Integer>> filter(ArrayList<HashMap> filters) {
        return this.filter(filters, null, null);
    }


    // -------------------------------------------------------------------------
    private boolean[] giniMask(Double minGini, Double minSpatialGini) {
        if (minGini != null) {
            if (minGini >= 1) {
                minGini = 1 - (1 - SarksUtilities.median(this.windGini)) * minGini;
            }
        }
        if (minSpatialGini != null) {
            if (minSpatialGini >= 1) {
                minSpatialGini =
                        1 - (1 - SarksUtilities.median(this.windGini)) * minSpatialGini;
            }
        }
        boolean[] mask = new boolean[this.windGini.length];
        for (int i=0; i<mask.length; i++) {mask[i] = false;}
        if (minGini != null) {
            for (int i=0; i<mask.length; i++) {
                if (this.windGini[i] < minGini) {mask[i] = true;}
            }
        }
        if (minSpatialGini != null && this.spatGini != null) {
            for (int i=0; i<mask.length; i++) {
                if (this.spatGini[i] < minSpatialGini) {mask[i] = true;}
            }
        }
        return mask;
    }
    
    public ArrayList<Float[][]> permutationDistribution(
            Integer reps, ArrayList<HashMap> filters,
            Long seed, Integer[][] permutations) {
        if (permutations == null) {
            permutations = SarksUtilities.generatePermutations(
                    reps, this.scores.length, seed);
        }
        Integer[][][] partPerm = SarksUtilities.partitionPermutations(
                permutations, this.nThreads);
        int hw0 = this.halfWindow;
        int sl0 = this.spatialLength;
        this.windowed = null;
        this.spatialWindowed = null;
        ArrayList<Float[][]> permResults = new ArrayList<Float[][]>();
        for (int r=0; r<permutations.length; r++) {
            permResults.add(new Float[filters.size()][2]);
        }
        for (int f=0; f<filters.size(); f++) {
            HashMap filt = filters.get(f);
            Integer hw = (Integer)filt.get("halfWindow");
            Double ming = (Double)filt.get("minGini");
            Integer sl = (Integer)filt.get("spatialLength");
            Double minsg = (Double)filt.get("minSpatialGini");
            if (hw != null && hw != this.halfWindow) {
                this.reset(hw, sl, false);
            } else if (sl != null && sl != this.spatialLength) {
                this.resetSpatial(sl, false);
            }
            boolean[] fmask = this.giniMask(ming, minsg);
            ArrayList<Future<Float[][]>> futures = new ArrayList<Future<Float[][]>>();
            ExecutorService pool = Executors.newFixedThreadPool(this.nThreads);
            for (int t=0; t<partPerm.length; t++) {
                Callable<Float[][]> tsarks =
                        new SarksPermuter(this, partPerm[t], fmask);
                Future<Float[][]> future = pool.submit(tsarks);
                futures.add(future);
            }
            ArrayList<Float[][]> unmerged = new ArrayList<Float[][]>();
            for (int t=0; t<partPerm.length; t++) {
                try {
                    unmerged.add(futures.get(t).get());
                } catch (Exception e) {}
            }
            pool.shutdown();
            Float[][] merged = SarksUtilities.mergeMaxResults(unmerged);
            for (int r=0; r<permutations.length; r++) {
                permResults.get(r)[f][0] = merged[r][0];
                permResults.get(r)[f][1] = merged[r][1];
            }
        }
        this.reset(hw0, sl0, true);
        return permResults;
    }
    public ArrayList<Float[][]> permutationDistribution(
            Integer reps, ArrayList<HashMap> filters, Long seed) {
        return this.permutationDistribution(reps, filters, seed, null);
    }
    public ArrayList<Float[][]> permutationDistribution(
            Integer reps, ArrayList<HashMap> filters, int seed) {
        return this.permutationDistribution(reps, filters, (long)seed, null);
    }
    public ArrayList<Float[][]> permutationDistribution(
            Integer reps, ArrayList<HashMap> filters) {
        return this.permutationDistribution(reps, filters, null, null);
    }


    public static class SarksPermuter implements Callable<Float[][]> {
        private Sarks base;
        private Integer[][] permutations;
        private boolean[] mask;
        public SarksPermuter(Sarks base0,
                             Integer[][] permutations0,
                             boolean[] mask0) {
            this.base = base0;
            this.permutations = permutations0;
            this.mask = mask0;
        }
        @Override
        public Float[][] call() {
            Float[][] out = new Float[permutations.length][2];
            for (int r=0; r<this.permutations.length; r++) {
                Sarks rsarks = new Sarks(this.base, permutations[r]);
                out[r][0] = SarksUtilities.max(rsarks.windowed, this.mask);
                out[r][1] = SarksUtilities.max(rsarks.spatialWindowed, this.mask);
            }
            return out;
        }
    }

    
    // -------------------------------------------------------------------------
    public static String printPermDist(ArrayList<HashMap> filters,
                                       ArrayList<Float[][]> permDist,
                                       String fileRoot,
                                       boolean spatial) throws Exception {
        String header = "halfWindow\tminGini\tspatialLength\tminSpatialGini\trep\tmax\n";
        BufferedWriter permWriter = null;
        StringBuilder sb = null;
        if (fileRoot != null) {
            if (!fileRoot.endsWith("/")) {fileRoot += "_";}
            String permFile = fileRoot + "windowed.tsv";
            permWriter = new BufferedWriter(new FileWriter(permFile));
            permWriter.write(header);
        } else {
            sb = new StringBuilder();
            sb.append(header);
        }
        int w = (spatial ? 1 : 0);
        for (int f=0; f<permDist.get(0).length; f++) {
            HashMap filt = filters.get(f);
            String hw = "" + (Integer)filt.get("halfWindow");
            String ming = "" + (Double)filt.get("minGini");
            String sl = "" + (Integer)filt.get("spatialLength");
            String minsg = "" + (Double)filt.get("minSpatialGini");
            String prefix = hw + "\t" + ming + "\t" + sl + "\t" + minsg + "\t";
            for (int r=0; r<permDist.size(); r++) {
                String line = prefix + r + "\t" + permDist.get(r)[f][w] + "\n";
                if (fileRoot != null) {
                    permWriter.write(line);
                } else {
                    sb.append(line);
                }
            }
        }
        if (fileRoot != null) {
            permWriter.close();
        } else {
            return sb.toString();
        }
        return null;
    }
    public static void printPermDists(ArrayList<HashMap> filters,
                                      ArrayList<Float[][]> permDist,
                                      String fileRoot) throws Exception {
        Sarks.printPermDist(filters, permDist, fileRoot, false);
        Sarks.printPermDist(filters, permDist, fileRoot+"_spatial", true);
    }

    public static String printThresholds(ArrayList<HashMap> filters,
                                         Float[][] thresholds,
                                         String fileRoot) throws Exception {
        String header = "halfWindow\tminGini\tspatialLength\tminSpatialGini\t";
        header += "theta\tspatialTheta\n";
        BufferedWriter thetaWriter = null;
        StringBuilder sb = null;
        if (fileRoot != null) {
            if (!fileRoot.endsWith("/")) {fileRoot += "_";}
            String thetaFile = fileRoot + "theta.tsv";
            thetaWriter = new BufferedWriter(new FileWriter(thetaFile));
            thetaWriter.write(header);
        } else {
            sb = new StringBuilder();
            sb.append(header);
        }
        for (int f=0; f<thresholds.length; f++) {
            HashMap filt = filters.get(f);
            String hw = "" + (Integer)filt.get("halfWindow");
            String ming = "" + (Double)filt.get("minGini");
            String sl = "" + (Integer)filt.get("spatialLength");
            String minsg = "" + (Double)filt.get("minSpatialGini");
            String prefix = hw + "\t" + ming + "\t" + sl + "\t" + minsg + "\t";
            String line = prefix + thresholds[f][0] + "\t" + thresholds[f][1] + "\n";
            if (fileRoot != null) {
                thetaWriter.write(line);
            } else {
                sb.append(line);
            }
        }
        if (fileRoot != null) {
            thetaWriter.close();
        } else {
            return sb.toString();
        }
        return null;
    }

    public String printPeaks(ArrayList<HashMap> filters,
                             Float[][] thresholds,
                             ArrayList<ArrayList<Integer>> iFilt,
                             String fileRoot,
                             int kmax) throws Exception {
        String header = "i\ts\tkmer\tkhat\tblock\twi\tgini\tscore\t" +
                        "windowed\tspatialWindowed\tkmax\t" +
                        "halfWindow\tminGini\ttheta\t" +
                        "spatialLength\tminSpatialGini\tspatialTheta\n";
        BufferedWriter peakWriter = null;
        StringBuilder sb = null;
        if (fileRoot != null) {
            if (!fileRoot.endsWith("/")) {fileRoot += "_";}            
            String peakFile = fileRoot + "peaks.tsv";
            peakWriter = new BufferedWriter(new FileWriter(peakFile));
            peakWriter.write(header);
        } else {
            sb = new StringBuilder();
            sb.append(header);
        }
        int hw0 = this.halfWindow;
        int sl0 = this.spatialLength;
        for (int f=0; f<thresholds.length; f++) {
            HashMap filt = filters.get(f);
            Integer hw = (Integer)filt.get("halfWindow");
            String ming = "" + (Double)filt.get("minGini");
            Integer sl = (Integer)filt.get("spatialLength");
            String minsg = "" + (Double)filt.get("minSpatialGini");
            String suffix = kmax + "\t" +
                            hw + "\t" + ming + "\t" + thresholds[f][0] + "\t" +
                            sl + "\t" + minsg + "\t" + thresholds[f][1] + "\n";
            if (hw != null && hw != this.halfWindow) {
                this.reset(hw, sl, true);
            } else if (sl != null && sl != this.spatialLength) {
                this.resetSpatial(sl, true);
            }
            for (int p=0; p<iFilt.get(f).size(); p++) {
                int i = iFilt.get(f).get(p);
                int s = this.sa[i];
                String km = this.kmer(s, kmax);
                float khat = (float)this.prefixAgreeSum(i, kmax);
                km = km.substring(0, java.lang.Math.round(khat));
                int blockIndex = this.sourceBlock(s);
                String block = this.transcripts[blockIndex];
                int wi = s - this.bounds[blockIndex];
                float gini = this.windGini[i];
                float score = (float)this.scores[blockIndex];
                Float sw = (this.spatialWindowed == null ? null : this.spatialWindowed[i]);
                String line = i + "\t" + s + "\t" + km + "\t" + khat + "\t" +
                              block + "\t" + wi + "\t" + gini + "\t" +
                              score + "\t" + this.windowed[i] + "\t" +
                              sw + "\t" + suffix;
                if (fileRoot != null) {
                    peakWriter.write(line);
                } else {
                    sb.append(line);
                }
            }
        }
        this.reset(hw0, sl0, true);
        if (fileRoot != null) {
            peakWriter.close();
        } else {
            return sb.toString();
        }
        return null;
    }
    public String printPeaks(ArrayList<HashMap> filters,
                             Float[][] thresholds,
                             ArrayList<ArrayList<Integer>> iFilt,
                             int kmax) throws Exception {
        return this.printPeaks(filters, thresholds, iFilt, null, kmax);
    }

    public String printBlockInfo(ArrayList<HashMap> filters,
                                 Float[][] thresholds,
                                 String block,
                                 String fileRoot,
                                 int kmax) throws Exception {
        ArrayList<Integer> eyes = new ArrayList<Integer>();
        int blockPos = this.transcriptPosition.get(block);
        for (int s=this.bounds[blockPos]; s<this.bounds[blockPos+1]; s++) {
            eyes.add(this.saInv[s]);
        }
        ArrayList<ArrayList<Integer>> iFilt =
                new ArrayList<ArrayList<Integer>>();
        for (int f=0; f<filters.size(); f++) {iFilt.add(eyes);}
        return this.printPeaks(filters, thresholds, iFilt, fileRoot, kmax);
    }


    // -------------------------------------------------------------------------
    public ArrayList<Integer> spatialSubPeaks(ArrayList<Integer> iFilt,
                                              Float theta,
                                              Double minGini) {
        if (this.spatialLength <= 1) {
            ArrayList<Integer> out = new ArrayList<Integer>();
            for (int i : iFilt) {out.add(this.sa[i]);}
            return out;
        }
        TreeSet<Integer> subPeaks = new TreeSet<Integer>();
        for (int i : iFilt) {
            int sLeft = this.sa[i];
            for (int s=sLeft; s<(sLeft+this.spatialLength); s++) {subPeaks.add(s);}
        }
        HashSet<Integer> filteredOut = new HashSet<Integer>();
        for (int s : subPeaks) {
            if (theta != null && this.windowed[this.saInv[s]] < theta) {
                filteredOut.add(s);
            } else if ((!subPeaks.contains(s-1)) &&
                       this.windowed[this.saInv[s-1]] >= theta) {
                filteredOut.add(s);
            }
        }
        subPeaks.removeAll(filteredOut);
        return new ArrayList<Integer>(subPeaks);
    }

    public ArrayList<int[]> mergeKmerIntervals(ArrayList<Integer> subpeaks, int kmax) {
        Collections.sort(subpeaks);
        int left = -1;
        int right = -1;
        int sLast = -2;
        ArrayList<int[]> out = new ArrayList<int[]>();
        for (int p=0; p<subpeaks.size(); p++) {
            int s = subpeaks.get(p);
            int khat = java.lang.Math.round(
                    (float)this.prefixAgreeSum(this.saInv[s], kmax));
            if ((s-1) == sLast) {
                right = java.lang.Math.max(right, s+khat);
            } else {
                if (left >= 0) {
                    out.add(new int[] {left, right});
                }
                left = s;
                right = s + khat;
                if (right > this.catSeq.length()) {
                    right = this.catSeq.length();
                }
            }
            sLast = s;
        }
        if (left >= 0) {
            out.add(new int[] {left, right});
        }
        return out;
    }

    public ArrayList<ArrayList<int[]>> multiMergeKmerIntervals(
            ArrayList<ArrayList<Integer>> iFilt,
            Float[][] thresholds,
            ArrayList<HashMap> filters,
            int kmax) {
        ArrayList<ArrayList<Integer>> subpeaks = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<int[]>> out = new ArrayList<ArrayList<int[]>>();
        int hw0 = this.halfWindow;
        int sl0 = this.spatialLength;
        for (int f=0; f<iFilt.size(); f++) {
            HashMap filt = filters.get(f);
            Integer hw = (Integer)filt.get("halfWindow");
            Double ming = (Double)filt.get("minGini");
            Integer sl = (Integer)filt.get("spatialLength");
            if ((hw != null) && (hw != this.halfWindow)) {
                this.reset(hw, sl, true);
            } else if ((sl != null) && (sl != this.spatialLength)) {
                this.resetSpatial(sl, true);
            }            
            subpeaks.add(this.spatialSubPeaks(iFilt.get(f),
                                              thresholds[f][1],
                                              ming));
            out.add(this.mergeKmerIntervals(subpeaks.get(f), kmax));
        }
        this.reset(hw0, sl0, true);
        return out;
    }

    public String[] mergedKmers(ArrayList<int[]> mergedKmerIntervals) {
        TreeSet<String> merged = new TreeSet<String>();
        for (int[] interval : mergedKmerIntervals) {
            merged.add(this.catSeq.substring(interval[0], interval[1]));
        }
        String[] out = new String[merged.size()];
        int kIndex = 0;
        for (String km : merged) {out[kIndex] = km; kIndex++;}
        return out;
    }
    public ArrayList<String[]> mergedKmers_(
            ArrayList<ArrayList<int[]>> mergedKmerIntervals) {
        ArrayList<String[]> out = new ArrayList<String[]>();
        for (int f=0; f<mergedKmerIntervals.size(); f++) {
            out.add(this.mergedKmers(mergedKmerIntervals.get(f)));
        }
        return out;
    }

    public String printMergedSubPeaks(ArrayList<HashMap> filters,
                                      Float[][] thresholds,
                                      ArrayList<ArrayList<int[]>> mergedKmerIntervals,
                                      String fileRoot,
                                      int kmax) throws Exception {
        String header = "i\ts\tkmer\tkhat\tblock\twi\tgini\tscore\t" +
                        "windowed\tspatialWindowed\tkmax\t" +
                        "halfWindow\tminGini\ttheta\t" +
                        "spatialLength\tminSpatialGini\tspatialTheta\n";
        BufferedWriter peakWriter = null;
        StringBuilder sb = null;
        if (fileRoot != null) {
            if (!fileRoot.endsWith("/")) {fileRoot += "_";}            
            String peakFile = fileRoot + "merged_peaks.tsv";
            peakWriter = new BufferedWriter(new FileWriter(peakFile));
            peakWriter.write(header);
        } else {
            sb = new StringBuilder();
            sb.append(header);
        }
        int hw0 = this.halfWindow;
        int sl0 = this.spatialLength;
        for (int f=0; f<thresholds.length; f++) {
            HashMap filt = filters.get(f);
            Integer hw = (Integer)filt.get("halfWindow");
            String ming = "" + (Double)filt.get("minGini");
            Integer sl = (Integer)filt.get("spatialLength");
            String minsg = "" + (Double)filt.get("minSpatialGini");
            String suffix = kmax + "\t" +
                            hw + "\t" + ming + "\t" + thresholds[f][0] + "\t" +
                            sl + "\t" + minsg + "\t" + thresholds[f][1] + "\n";
            if (hw != null && hw != this.halfWindow) {
                this.reset(hw, sl, true);
            } else if (sl != null && sl != this.spatialLength) {
                this.resetSpatial(sl, true);
            }
            for (int p=0; p<mergedKmerIntervals.get(f).size(); p++) {
                int[] sInterval = mergedKmerIntervals.get(f).get(p);
                int i = this.saInv[sInterval[0]];
                String km = this.catSeq.substring(sInterval[0], sInterval[1]);
                float khat = (float)this.prefixAgreeSum(i, kmax);
                int blockIndex = this.sourceBlock(sInterval[0]);
                String block = this.transcripts[blockIndex];
                int wi = sInterval[0] - this.bounds[blockIndex];
                float gini = this.windGini[i];
                float score = (float)this.scores[blockIndex];
                Float sw = (this.spatialWindowed == null ? null : this.spatialWindowed[i]);
                String line = i + "\t" + sInterval[0] + "\t" +
                              km + "\t" + khat + "\t" +
                              block + "\t" + wi + "\t" + gini + "\t" +
                              score + "\t" + this.windowed[i] + "\t" +
                              sw + "\t" + suffix;
                if (fileRoot != null) {
                    peakWriter.write(line);
                } else {
                    sb.append(line);
                }
            }
        }
        this.reset(hw0, sl0, true);
        if (fileRoot != null) {
            peakWriter.close();
        } else {
            return sb.toString();
        }
        return null;
    }
    public String printMergedSubPeaks(ArrayList<HashMap> filters,
                                      Float[][] thresholds,
                                      ArrayList<ArrayList<int[]>> mergedKmerIntervals,
                                      int kmax) throws Exception {
        return this.printMergedSubPeaks(
                filters, thresholds, mergedKmerIntervals, null, kmax);
    }


    // -------------------------------------------------------------------------
    public int getHalfWindow() {return this.halfWindow;}
    public int getSpatialLength() {return this.spatialLength;}
    public String getCatSeq() {return this.catSeq;}
    public double[] getScores() {return this.scores;}
    public String[] getTranscripts() {return this.transcripts;}
    public HashMap<String,Integer> getTranscriptPosition() {return this.transcriptPosition;}
    public int[] getSuffixArray() {return this.sa;}
    public int[] getInverseSuffixArray() {return this.saInv;}
    public float[] getGini() {return this.windGini;}
    public float[] getSpatialGini() {return this.spatGini;}
    public float[] getYhat() {return this.windowed;}
    public float[] getYdoubleHat() {return this.spatialWindowed;}

    public int s2i(int s) {return this.saInv[s];}
    public int i2s(int i) {return this.sa[i];}
}
