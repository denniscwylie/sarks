package dcw.sarks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

public class SarksUtilities {

    public static Float max(float[] arr, boolean[] mask) {
        if ((arr == null) || (arr.length == 0)) {return null;}
        float out = Float.NEGATIVE_INFINITY;
        for (int i=0; i<arr.length; i++) {
            if ((!mask[i]) && (arr[i] > out)) {out = arr[i];}
        }
        return out;
    }

    public static int max(int[] arr, int len) {
        int out = Integer.MIN_VALUE;
        for (int i=0; i<len; i++) {if (arr[i] > out) {out = arr[i];}}
        return out;
    }

    public static Float mean(Float[] arr) {
        if ((arr == null) || (arr.length == 0)) {return null;}
        if (arr[0] == null) {return null;}
        double out = 0;
        for (int i=0; i<arr.length; i++) {out += arr[i];}
        return (float)out / (float)arr.length;
    }

    public static Float sd(Float[] arr) {
        if ((arr == null) || (arr.length == 0)) {return null;}
        if (arr[0] == null) {return null;}        
        double sum = 0;
        double sumsq = 0;
        for (int i=0; i<arr.length; i++) {
            sum += arr[i];
            sumsq += (arr[i] * arr[i]);
        }
        double mean = sum / (double)arr.length;
        double meansq = sumsq / (double)arr.length;
        return (float)java.lang.Math.sqrt(
            ((double)(arr.length) / (double)(arr.length-1)) *
            (meansq - (mean*mean))
        );
    }
        
    public static float median(float[] arr) {
        float[] arrToSort = new float[arr.length];
        for (int i=0; i<arr.length; i++) {arrToSort[i] = arr[i];}
        java.util.Arrays.sort(arrToSort);
        int mid1 = arr.length / 2;
        int mid2 = mid1;
        if (arr.length > 0 && arr.length % 2 == 0) {mid2--;}
        return (float)((arrToSort[mid1] + arrToSort[mid2]) / 2.0);
    }


    // -------------------------------------------------------------------------
    public static Integer[][] generatePermutations(int reps, int n, Long seed) {
        Random rand = new Random();
        if (seed != null) {rand.setSeed(seed);}
        Integer[][] permutations = new Integer[reps][n];
        for (int r=0; r<reps; r++) {
            ArrayList<Integer> shuffled = new ArrayList<Integer>();
            for (int b=0; b<n; b++) {shuffled.add(b);}
            Collections.shuffle(shuffled, rand);
            for (int b=0; b<n; b++) {permutations[r][b] = shuffled.get(b);}
        }
        return permutations;
    }

    public static Integer[][][] partitionPermutations(
            Integer[][] permutations, int partitions) {
        Integer[][][] partitioned = new Integer[partitions][][];
        int minSize = permutations.length / partitions;
        int nBigger = permutations.length % partitions;
        int rOffset = 0;
        for (int p=0; p<partitions; p++) {
            if (p < nBigger) {
                partitioned[p] = new Integer[minSize+1][];
            } else {
                partitioned[p] = new Integer[minSize][];
            }
            for (int r=0; r<partitioned[p].length; r++) {
                partitioned[p][r] = permutations[r+rOffset];
            }
            rOffset += partitioned[p].length;
        }
        return partitioned;
    }

    public static Float[][] mergeMaxResults(ArrayList<Float[][]> results) {
        int reps = 0;
        for (int t=0; t<results.size(); t++) {reps += results.get(t).length;}
        Float[][] out = new Float[reps][2];
        int rOffset = 0;
        for (int t=0; t<results.size(); t++) {
            Float[][] tres = results.get(t);
            for (int r=0; r<tres.length; r++) {
                out[r+rOffset][0] = tres[r][0];
                out[r+rOffset][1] = tres[r][1];
            }
            rOffset += tres.length;
        }
        return out;
    }

    public static Float[][] thresholdsFromPermutations(
            ArrayList<Float[][]> permDists, float nSigma) {
        Float[][] thresholds = new Float[permDists.get(0).length][2];
        for (int f=0; f<thresholds.length; f++) {
            Float[] winMaxes = new Float[permDists.size()];
            Float[] spatMaxes = new Float[permDists.size()];
            for (int r=0; r<winMaxes.length; r++) {
                winMaxes[r] = permDists.get(r)[f][0];
                spatMaxes[r] = permDists.get(r)[f][1];
            }
            if (winMaxes[0] != null) {
                thresholds[f][0] = SarksUtilities.mean(winMaxes) +
                                   (nSigma * SarksUtilities.sd(winMaxes));
            } else {
                thresholds[f][0] = null;
            }
            if (spatMaxes[0] != null) {
                thresholds[f][0] = null;
                thresholds[f][1] = SarksUtilities.mean(spatMaxes) +
                                   (nSigma * SarksUtilities.sd(spatMaxes));
            } else {
                thresholds[f][1] = null;
            }
        }
        return thresholds;
    }
    public static Float[][] thresholdsFromPermutations(
            ArrayList<Float[][]> permDists, double nSigma) {
        return SarksUtilities.thresholdsFromPermutations(permDists,
                                                         (float)nSigma);
    }

    public static int falsePositives(
            ArrayList<Float[][]> permDistsTest, Float[][] thresholds) {
        int[] fp = new int[permDistsTest.size()];
        for (int f=0; f<thresholds.length; f++) {
            for (int r=0; r<permDistsTest.size(); r++) {
                if (((thresholds[f][0] != null) &&
                     (permDistsTest.get(r)[f][0] >= thresholds[f][0])) ||
                    ((thresholds[f][1] != null) &&
                     (permDistsTest.get(r)[f][1] >= thresholds[f][1]))) {
                    fp[r] = 1;
                }
            }
        }
        int out = 0;
        for (int r=0; r<fp.length; r++) {out += fp[r];}
        return out;
    }

    
    // -------------------------------------------------------------------------
    public static HashMap<String,String> readFasta(String filename) throws Exception {
        HashMap<String,String> out = new HashMap<String,String>();
        File file = new File(filename);        
        if (file.exists()) {
            String curName = null;
            StringBuilder sb = new StringBuilder();
            BufferedReader in = null;
            if (filename.toLowerCase().endsWith(".gz") ||
                filename.toLowerCase().endsWith(".gzip")) {
                in = new BufferedReader(new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(filename))));
            } else {
                in = new BufferedReader(new FileReader(filename));
            }
            String curLine = in.readLine();
            while (curLine != null) {
                if (curLine.startsWith(">")) {
                    if (curName != null) {
                        out.put(curName, sb.toString());
                        sb = new StringBuilder();
                    }
                    curName = curLine.substring(1);
                } else {
                    sb.append(curLine);
                }
                curLine = in.readLine();
            }
            if (curName != null && !out.containsKey(curName)) {
                out.put(curName, sb.toString());
            }
        }
        return out;
    }

    public static TreeMap<String,Double> readScores(String filename) throws Exception {
        TreeMap<String,Double> out = new TreeMap<String,Double>();
        File file = new File(filename);
        if (file.exists()) {
            BufferedReader in = null;
            if (filename.toLowerCase().endsWith(".gz") ||
                filename.toLowerCase().endsWith(".gzip")) {
                in = new BufferedReader(new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(filename))));
            } else {
                in = new BufferedReader(new FileReader(filename));
            }
            String curLine = in.readLine();
            while (curLine != null) {
                String[] tokens = curLine.split("\t");
                try {
                    double score = Double.parseDouble(tokens[1]);
                    out.put(tokens[0], score);
                } catch (NumberFormatException nfe) {}
                curLine = in.readLine();                
            }
        }
        return out;
    }
}
