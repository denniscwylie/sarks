package dcw.sarks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class SarksMain {
    public static void sarkselect(String fasta,
                                  String scoreFile,
                                  int[] halfWindows,
                                  int[] spatialLengths,
                                  double[] minGinis,
                                  int reps,
                                  Integer kMax,
                                  double nSigma,
                                  String outDir,
                                  Integer nThreads,
                                  Long seed) throws Exception {
        if (!outDir.endsWith("/")) {outDir += "/";}
        if (spatialLengths == null) {spatialLengths = new int[] {0};}
        if (minGinis == null) {minGinis = new double[] {Double.NEGATIVE_INFINITY};}
        if (kMax == null) {kMax = 12;}
        if (nThreads == null) {nThreads = 1;}
        ArrayList<HashMap> filters = new ArrayList<HashMap>();
        for (int hw : halfWindows) {
            for (int sl : spatialLengths) {
                for (double mg : minGinis) {
                    HashMap pars = new HashMap();
                    pars.put("halfWindow", hw);
                    pars.put("minGini", mg);
                    pars.put("spatialLength", sl);
                    pars.put("minSpatialGini", mg);
                    filters.add(pars);
                }
            }
        }
        Sarks sarks = new Sarks(fasta,
                                scoreFile,
                                (int)filters.get(0).get("halfWindow"),
                                (int)filters.get(0).get("spatialLength"),
                                nThreads,
                                false);
        ArrayList<Float[][]> permDist =
                sarks.permutationDistribution(reps, filters, seed);
        (new File(outDir)).mkdirs();
        Sarks.printPermDists(filters, permDist, outDir+"permdists");
        Float[][] thresholds =
                SarksUtilities.thresholdsFromPermutations(permDist, nSigma);
        Sarks.printThresholds(filters, thresholds, outDir);
        ArrayList<ArrayList<Integer>> eyes = sarks.filter(filters, thresholds, true);
        sarks.printPeaks(filters, thresholds, eyes, outDir, kMax.intValue());
        ArrayList<ArrayList<int[]>> mergedKmerIntervals =
                sarks.multiMergeKmerIntervals(eyes, thresholds, filters, kMax);
        sarks.printMergedSubPeaks(filters, thresholds,
                                  mergedKmerIntervals, outDir, kMax);
    }

    public static ArrayList<HashMap> readThresholds(String thresholdFile) throws Exception {
        ArrayList<HashMap> filters = new ArrayList<HashMap>();
        BufferedReader in = new BufferedReader(new FileReader(thresholdFile));
        String curLine = in.readLine();
        curLine = in.readLine();   // ignore header line
        while (curLine != null) {
            HashMap pars = new HashMap();
            String[] entries = curLine.split("\t");
            pars.put("halfWindow", Integer.parseInt(entries[0]));
            pars.put("minGini", Double.parseDouble(entries[1]));
            pars.put("spatialLength", Integer.parseInt(entries[2]));
            if (entries[3].equals("null")) {
                pars.put("minSpatialGini", null);
            } else {
                pars.put("minSpatialGini", Double.parseDouble(entries[3]));
            }
            if (entries[4].equals("null")) {
                pars.put("theta", null);
            } else {
                pars.put("theta", Float.parseFloat(entries[4]));
            }
            if (entries[5].equals("null")) {
                pars.put("spatialTheta", null);
            } else {
                pars.put("spatialTheta", Float.parseFloat(entries[5]));
            }
            filters.add(pars);
            curLine = in.readLine();
        }
        in.close();
        return filters;
    }

    public static int sarkstest(String fasta,
                                String scoreFile,
                                String thresholdFile,
                                int testReps,
                                String outDir,
                                Integer nThreads,
                                Long seed) throws Exception {
        if (!outDir.endsWith("/")) {outDir += "/";}
        if (nThreads == null) {nThreads = 1;}
        ArrayList<HashMap> filters = SarksMain.readThresholds(thresholdFile);
        Float[][] thresholds = new Float[filters.size()][2];
        for (int f=0; f<thresholds.length; f++) {
            thresholds[f][0] = (Float)filters.get(f).get("theta");
            thresholds[f][1] = (Float)filters.get(f).get("spatialTheta");
        }
        Sarks sarks = new Sarks(fasta,
                                scoreFile,
                                (int)filters.get(0).get("halfWindow"),
                                (int)filters.get(0).get("spatialLength"),
                                nThreads,
                                false);
        ArrayList<Float[][]> permDistTest =
                sarks.permutationDistribution(testReps, filters, seed);
        (new File(outDir)).mkdirs();
        Sarks.printPermDists(filters, permDistTest, outDir+"test_permdists");
        return SarksUtilities.falsePositives(permDistTest, thresholds);
    }

    public static String helpMessage(String application) {
        if (application == null) {
            return "sarks (Suffix Array Kernel Smoothing)\n" +
                "  https://academic.oup.com/bioinformatics/article/35/20/3944/5418797\n" +
                "  https://github.com/denniscwylie/sarks\n" +
                "\nusage: java [jvm options] -jar sarks.jar <command> [arguments]\n" +
                "\nCommands:\n" +
                "  select    main step for motif discovery\n" +
                "  test      estimate false positive rate associated with\n" +
                "            motif set determined by select command\n";
        } else if (application.equals("select")) {
            return "usage:\njava [jvm options] -jar sarks.jar select [-h] [-f FASTA] [-s SCORES]\n" +
                "                                         [-w HALFWINDOW] [-l SPATIALLENGTH]\n" +
                "                                         [-g MINGINI] [-r REPS] [-z NSIGMA]\n" +
                "                                         [-e SEED] [-t NTHREADS] [-o OUTDIR]\n" +
                "\nrequired arguments:\n" +
                "  -f, --fasta FILE        input fasta file containing sequences to analyze\n" +
                "                          (can be gzipped)\n" +
                "  -s, --scores FILE       input scores tsv file:\n" +
                "                          col1=seqids, col2=numeric scores\n" +
                "                          (can be gzipped)\n" +
                "  -w, --halfwindow INT    half window width (kappa) for first smoothing pass,\n" +
                "                          can supply multiple values using commas (no spaces)\n" +
                "  -o, --outdir DIR        output directory be created/overwritten\n" +
                "\noptional arguments:\n" +
                "  -l, --spatiallength INT spatial smoothing length (lambda), can supply\n" +
                "                          multiple values using commas (no spaces)\n" +
                "                          [default 0]\n" +
                "  -g, --mingini FLOAT     parameter for calculation of Gini filter (gamma),\n" +
                "                          can supply multiple values using commas (no spaces)\n" +
                "                          [default 1.1]\n" +
                "  -r, --reps INT          number R of permutations used to set significance\n" +
                "                          thresholds  [default 100]\n" +
                "  -z, --nsigma FLOAT      multiple z of standard deviations above mean (of\n" +
                "                          max smoothed suffix scores obtained after randomly\n" +
                "                          permuting sequence scores) defining threshold\n" +
                "                          (Section S2.6, Eq (S24-S25) of paper)\n" +
                "                          [default 4.0]\n" +
                "  -e, --seed INT          seed for random number generator\n" +
                "  -t, --threads INT       number of threads to use for permutational analyses\n" +
                "                          [default 1]\n" +
                "  -h, --help              show this help message and exit\n";
        } else if (application.equals("test")) {
            return "usage:\njava [jvm options] -jar sarks.jar test [-h] [-f FASTA] [-s SCORES]\n" +
                "                                       [-i INDIR] [-r REPS]\n" +
                "                                       [-e SEED] [-t NTHREADS]\n" +
                "\nrequired arguments:\n" +
                "  -f, --fasta FILE        input fasta file containing sequences to analyze\n" +
                "                          (can be gzipped)\n" +
                "  -s, --scores FILE       input scores tsv file:\n" +
                "                          col1=seqids, col2=numeric scores\n" +
                "                          (can be gzipped)\n" +
                "  -i, --indir DIR         input directory for test command is output directory\n" +
                "                          from select command\n" +
                "\noptional arguments:\n" +
                "  -r, --reps INT          number R_2 of permutations used to test\n" +
                "                          significance thresholds  [default 100]\n" +
                "  -e, --seed INT          seed for random number generator;\n" +
                "                          should NOT be the same seed as was used\n" +
                "                          in select command used to generate input directory\n" +
                "  -t, --threads INT       number of threads to use during permutational analyses\n" +
                "                          [default 1]\n" +
                "  -h, --help              show this help message and exit\n";
        }
        return null;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.out.println(SarksMain.helpMessage(null));
            System.exit(0);
        }
        String application = args[0].toLowerCase();
        String key = null;
        if (application.equalsIgnoreCase("-h") ||
            application.equalsIgnoreCase("--help")) {
            System.out.println(SarksMain.helpMessage(null));
            System.exit(0);
        }
        if (args.length <= 2) {
            if ((args.length == 1) ||
                args[1].equalsIgnoreCase("-h") ||
                args[1].equalsIgnoreCase("--help")) {
                System.out.println(SarksMain.helpMessage(application));
                System.exit(0);
            }
        }
        HashMap options = new HashMap();
        // set defaults
        options.put("spatialLength", (Integer)0);
        options.put("minGini", (Float)1.1f);
        options.put("reps", (Integer)100);
        options.put("nSigma", (Double)4.0);
        // end set defaults
        for (String arg : args) {
            if (arg.startsWith("-")) {
                key = arg;
                while (key.startsWith("-")) {key = key.substring(1);}
            } else if (key != null) {
                String lkey = key.toLowerCase();
                if (lkey.equals("w") || lkey.equals("halfwindow") ||
                    lkey.equals("l") || lkey.equals("spatiallength")) {
                    String[] tokens = arg.split(",");
                    int[] values = new int[tokens.length];
                    for (int t=0; t<tokens.length; t++) {
                        values[t] = Integer.parseInt(tokens[t]);
                    }
                    if (lkey.equals("w") || lkey.equals("halfwindow")) {
                        options.put("halfWindow", values);
                    } else if (lkey.equals("l") || lkey.equals("spatiallength")) {
                        options.put("spatialLength", values);
                    }
                } else if (lkey.equals("g") || lkey.equals("mingini")) {
                    String[] tokens = arg.split(",");
                    double[] values = new double[tokens.length];
                    for (int t=0; t<tokens.length; t++) {
                        values[t] = Double.parseDouble(tokens[t]);
                    }
                    options.put("minGini", values);
                } else if (lkey.equals("i") || lkey.equals("indir")) {
                    options.put("inDir", arg);
                } else if (lkey.equals("f") || lkey.equals("fasta")) {
                    options.put("fasta", arg);
                } else if (lkey.equals("s") || lkey.equals("scores")) {
                    options.put("scoreFile", arg);
                } else if (lkey.equals("r") || lkey.equals("reps")) {
                    options.put("reps", Integer.parseInt(arg));
                } else if (lkey.equals("o") || lkey.equals("outdir")) {
                    options.put("outDir", arg);
                } else if (lkey.equals("z") || lkey.equals("nsigma")) {
                    options.put("nSigma", Double.parseDouble(arg));
                } else if (lkey.equals("e") || lkey.equals("seed")) {
                    options.put("seed", Long.parseLong(arg));
                } else if (lkey.equals("k") || lkey.equals("kmax")) {
                    options.put("kMax", Integer.parseInt(arg));
                } else if (lkey.equals("t") ||
                           lkey.equals("threads") ||
                           lkey.equals("nthreads")) {
                    options.put("nThreads", Integer.parseInt(arg));
                }
            }
        }
        if ((options.get("fasta") == null) ||
            (options.get("scoreFile") == null)) {
            System.out.println("Must specify -f FASTA and -s SCORES input files.");
            System.exit(-1);
        } else {
            if (!(new File((String)options.get("fasta"))).exists()) {
                System.out.println((String)options.get("fasta") + " not found.");
                System.exit(-1);
            }
            if (!(new File((String)options.get("scoreFile"))).exists()) {
                System.out.println((String)options.get("scoreFile") + " not found.");
                System.exit(-1);
            }
        }
        if (application.equals("select")) {
            SarksMain.sarkselect((String)options.get("fasta"),
                                 (String)options.get("scoreFile"),
                                 (int[])options.get("halfWindow"),
                                 (int[])options.get("spatialLength"),
                                 (double[])options.get("minGini"),
                                 (int)options.get("reps"),
                                 (Integer)options.get("kMax"),
                                 (double)options.get("nSigma"),
                                 (String)options.get("outDir"),
                                 (Integer)options.get("nThreads"),
                                 (Long)options.get("seed"));
        } else if (application.equals("test")) {
            String thresholdFile = (String)options.get("inDir");
            if (!thresholdFile.endsWith("theta.tsv")) {
                if (!thresholdFile.endsWith("/")) {
                    thresholdFile += "/";
                }
                if (options.get("outDir") == null) {
                    options.put("outDir", thresholdFile);
                }
                thresholdFile += "theta.tsv";
            }
            int fp = SarksMain.sarkstest((String)options.get("fasta"),
                                         (String)options.get("scoreFile"),
                                         thresholdFile,
                                         (int)options.get("reps"),
                                         (String)options.get("outDir"),
                                         (Integer)options.get("nThreads"),
                                         (Long)options.get("seed"));
            System.out.println(""+fp);
        }
    }
}
