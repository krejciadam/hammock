/**
 * Main class. Parses command line arguments and runs appropriate routines
 */
package hammock;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 *
 * @author Adam Krejci
 */
public class Hammock {

    public static final char separatorChar = File.separatorChar; //system separator
    private static final String parentDir = (new File(Hammock.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getParentFile().getPath());

    //common
    private static String inputFileName = null;
    public static String workingDirectory = null;
    public static int nThreads = 4;
    private static List<String> labels = null;

    public static List<String> getLabels() {
        return labels;
    }
    public static double minCorrelation = -1.0;
    private static String labelString = null;
    private static Boolean fullClustering = null;
    public static ExecutorService threadPool = null;
    public static boolean inGalaxy = false;

    //common - preset
    private static String initialAlignedCl;
    private static String initialAlignedSequenceCsv;
    private static String initialAlignedClusterCsv;
    private static String finalCl;
    private static String finalRemainingSequences;
    private static String finalSequenceCsv;
    private static String finalClustersCsv;
    private static String initialAlnFolder;
    private static String inputStatistics;
    private static String finalAlnFolder;
    public static Logger logger;

    //greedy
    private static String inputType = "fasta";
    private static String matrixFile = parentDir + separatorChar + "matrices" + separatorChar + "blosum62.txt";
    private static Integer greedyThreshold = null;
    private static int shiftPenalty = 0;
    private static int maxShift = 3;
    private static int[] ignore = new int[]{1};

    public static Map<Integer, UniqueSequence> pivotMap = null;                 //this should be solved another way, without global variables
    public static Map<Integer, List<AligningScorerResult>> foundSeqsMap = null; //this should be solved another way, without global variables
    //hmm
    private static Integer sizeThreshold = null;
    private static Integer countThreshold = null;
    private static Double partThreshold = null;
    private static String assignThresholdSequenceString = null;
    private static String overlapThresholdSequenceString = null;
    private static String hhMergeThresholdSequenceString = null;
    private static double[] assignThresholdSequence = null;
    private static double[] overlapThresholdSequence = null;
    private static double[] hhMergeThresholdSequence = null;
    private static boolean[] fullHHClustering = null;
    public static boolean relativeHHScore = false;
    public static boolean relativeHmmScore = false;
    public static Integer minMatchStates = 4;
    public static Double minIc = 1.2;
    public static Double maxGapProportion = 0.05;
    public static Integer maxAlnLength = null;
    public static int maxInnerGaps = 0;
    public static boolean extensionIncreaseLength = false;

    //galaxy
    private static String galaxyFinalClustersCsv;
    private static String galaxyFinalSequenceCsv;

    public static void main(String[] args) throws IOException, HammockException, InterruptedException, ExecutionException, Exception {
        try {
            parseArgs(args);
        } catch (FileFormatException e) {
            logger.logAndStderr("Error. Probably wrong input file format? Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            logger.logAndStderr(errors.toString());
            if (inGalaxy) {
                System.err.println("Error. Probably wrong input file format. Trace: \n");
                System.err.println(errors.toString());
            }
        } catch (NullPointerException e) {
            logger.logAndStderr("Error. Maybe wrong input file format? Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            logger.logAndStderr(errors.toString());
            if (inGalaxy) {
                System.err.println("Error. Maybe wrong input file format? Trace: \n");
                System.err.println(errors.toString());
            }
        } catch (DataException e) {
            logger.logAndStderr("Error. Maybe wrong input file format or wrong set of labels? Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            logger.logAndStderr(errors.toString());
            if (inGalaxy) {
                System.err.println("Error. Maybe wrong input file format or wrong set of labels? Trace: \n");
                System.err.println(errors.toString());
            }
        } catch (Exception e) {
            logger.logAndStderr("Error. Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            logger.logAndStderr(errors.toString());
            if (inGalaxy) {
                System.err.println("Error. Trace: \n");
                System.err.println(errors.toString());
            }
        }
        if (threadPool != null) {
            threadPool.shutdown();
        }
    }

    /**
     * Parses command line arguments and runs appropriate routines
     *
     * @param args
     * @throws IOException
     * @throws FileFormatException
     * @throws HammockException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws Exception
     */
    private static void parseArgs(String[] args) throws IOException, FileFormatException, HammockException, InterruptedException, ExecutionException, Exception {
        if (args.length < 1) {
            System.err.println("First argument must be a command. Run \"java -jar Hammpck.jar --help\" for more info.");
            return;
        }

        if (args[0].equals("--help") || args[0].equals("help") || args[0].equals("-h")) {
            displayHelp();
            return;
        }

        StringBuilder argsString = new StringBuilder();
        for (String arg : args) {
            argsString.append(" ").append(arg);
        }

        if (args[0].equals("cluster")) {
            fullClustering = false;
            parseCommonArgs(args);
            if (inGalaxy) {
                throw new DataException("Error. Can't run in cluster mode in Galaxy");
            }
            parseClusteringArgs(args);
            if ((checkCommonArgs()) && (checkClusteringArgs())) {
                logger.logWithTime("Program started in mode \"cluster\".");
                logger.logWithoutTime("Parameters: " + "\n" + argsString.toString() + "\n");
                runClustering();
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }

        if (args[0].equals("greedy")) {
            fullClustering = false;
            parseCommonArgs(args);
            if (inGalaxy) {
                throw new DataException("Error. Can't run in greedy mode in Galaxy");
            }
            parseGreedyArgs(args);
            if (checkCommonArgs() && checkGreedyArgs()) {
                logger.logWithTime("Program started in mode \"greedy\".");
                logger.logWithoutTime("Parameters: " + "\n" + argsString.toString() + "\n");
                runGreedyClustering();
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }

        if (args[0].equals("full")) {
            fullClustering = true;
            parseCommonArgs(args);
            parseGreedyArgs(args);
            parseClusteringArgs(args);
            if (checkCommonArgs() && checkGreedyArgs() && checkClusteringArgs()) {
                logger.logWithTime("Program started in mode \"full\".");
                logger.logWithoutTime("Parameters: " + "\n" + argsString.toString() + "\n");
                runFull();
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }

        System.err.println("First argument must be a command. Run \"java -jar Hammpck.jar --help\" for more info.");
    }

    /**
     * Displays help in command line
     */
    private static void displayHelp() {
        System.err.println("HAMMOCK - a Hidden Markov Model based peptide sequence clustering tool. For more details and examples, see manual at : http://www.recamo.cz/en/software/hammock-cluster-peptides/");
        System.err.println("\nVersion: 1.0.0\n");
        System.err.println("\n------MANUAL------\n");
        System.err.println("Synopsis: java -jar Hammock.jar <mode> <param1> <param2> <param3>...");
        System.err.println("\nModes: \n");
        System.err.println("full\ngreedy\ncluster");
        System.err.println("\n------parameters common for all modes------\n");
        System.err.println("-i, --input <file>   path to input file");
        System.err.println("-d, --outputDirectory <directory>   directory for all output files");
        System.err.println("-t, --threads <int>   number of threads to be used");
        System.err.println("-l, --labels <str,str,str...>  list of sequence labels to use");
        System.err.println("\n------parameters specific for greedy mode------\n");
        System.err.println("-f, --file_format <[fasta,seq,tab]>   file format of input file specified by -i");
        System.err.println("-m, --matrix <file>   path to substitution matrix file");
        System.err.println("-g, --greedy_threshold 〈int〉   Minimal sequence neede for a sequence to join a cluster");
        System.err.println("-x, --max_shift 〈int〉   maximal sequence-sequence shift. Nonnegative int");
        System.err.println("-p, --gap_penalty 〈int〉   penalty for each position of shift. Nonpositive int");
        System.err.println("\n------parameters specific for cluster mode------\n");
        System.err.println("-a, --part_threshold <float [0, 1]>   the proportion of clusters to be used as cores");
        System.err.println("-s, --size_threshold 〈int〉   minimal size of a cluster to be used a core");
        System.err.println("-c, --count_threshold 〈int〉   the exact number of clusters to be used as cores");
        System.err.println("-n, --assign_thresholds 〈float,float,float...〉   sequence of threshold values for sequence to cluster assignment");
    }

    /**
     * Runs full routine - greedy clustering followed by hmm clustering
     *
     * @throws IOException
     * @throws HammockException
     * @throws FileFormatException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws Exception
     */
    private static void runFull() throws IOException, HammockException, FileFormatException, InterruptedException, ExecutionException, Exception {
        if (inGalaxy) {
            initialAlignedCl = Settings.getInstance().getTempDirectory() + separatorChar + "greedy_aligned_cl";
            finalClustersCsv = galaxyFinalClustersCsv;
            finalSequenceCsv = galaxyFinalSequenceCsv;
        }
        runGreedyClustering();
        inputFileName = initialAlignedCl;
        runClustering();
    }

    /**
     * Runs greedy clustering routine
     *
     * @throws IOException
     * @throws FileFormatException
     * @throws HammockException
     * @throws InterruptedException
     * @throws ExecutionException
     */
    private static void runGreedyClustering() throws IOException, FileFormatException, HammockException, InterruptedException, ExecutionException {

        List<UniqueSequence> sequences = null;
        logger.logAndStderr("Loading input sequences...");
        if (inputType.equals("fasta")) {
            sequences = FileIOManager.loadUniqueSequencesFromFasta(inputFileName);
        }
        if (inputType.equals("seq")) {
            sequences = FileIOManager.loadUniqueSequencesFromFile(inputFileName);
        }
        if (inputType.equals("tab")) {
            sequences = FileIOManager.loadUniqueSequencesFromTable(inputFileName, ignore);
        }

        if (sequences == null) {
            throw new HammockException("Error, this should have been checked.");
        }
        logger.logAndStderr(sequences.size() + " unique sequences loaded.");
        logger.logAndStderr(Hammock.getTotalSequenceCount(sequences) + " total sequences loaded.");

        if (labels != null) {
            sequences = filterSequencesForLabels(sequences, labels);
        }

        if (labels == null) {
            labels = getSortedLabels(sequences);
        }
        
        logger.logAndStderr(sequences.size() + " unique sequences after non-specified labels filtered out");
        logger.logAndStderr(new Cluster(sequences, -1).size() + " total sequences after non-specified labels fileterd out");

        if (sequences.isEmpty()) {
            throw new FileFormatException("Error. No sequences (with specified labels) to cluster.");
        }
        if (greedyThreshold == null) {
            greedyThreshold = setGreedyThreshold(sequences);
            logger.logAndStderr("Greedy clustering threshold not set. Setting automatically based on mean sequence length to: " + greedyThreshold);
        }
        logger.logAndStderr("Generating input statistics...");
        if (!inGalaxy) {
            FileIOManager.saveInputStatistics(sequences, labels, inputStatistics);
        }
        SequenceClusterer clusterer;
        ShiftedScorer scorer = new ShiftedScorer(FileIOManager.loadScoringMatrix(matrixFile), shiftPenalty, maxShift);

        clusterer = new AligningGreedySequenceClusterer(greedyThreshold);

        logger.logAndStderr("Greedy clustering...");
        Long time = System.currentTimeMillis();
        List<Cluster> clusters = clusterer.cluster(sequences, scorer);
        logger.logAndStderr("Ready. Clustering time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Resulting clusers: " + clusters.size());
        logger.logAndStderr("Building MSAs...");
        for (Cluster cl : clusters) {
            FileIOManager.makeShiftedClusterAlignment(pivotMap.get(cl.getId()), foundSeqsMap.get(cl.getId()), cl.getId());
            cl.setAsHasMSA();
        }
        logger.logAndStderr("Ready. Total time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Saving results to output files...");
        FileIOManager.saveClustersToFile(clusters, initialAlignedCl);
        if (!inGalaxy) {
            FileIOManager.saveClusterSequencesToCsv(clusters, initialAlignedSequenceCsv, labels);
            FileIOManager.SaveClustersToCsv(clusters, initialAlignedClusterCsv, labels);
        }
        logger.logAndStderr("Greedy clustering results in: " + initialAlignedSequenceCsv);
        logger.logAndStderr("and: " + initialAlignedClusterCsv);
        logger.logAndStderr("and: " + initialAlignedCl);
    }

    /**
     * Runs hmm clustering routine
     *
     * @throws IOException
     * @throws FileFormatException
     * @throws ExecutionException
     * @throws Exception
     */
    private static void runClustering() throws IOException, FileFormatException, ExecutionException, Exception {
        logger.logAndStderr("\nLoading clusters...");
        List<Cluster> inClusters = FileIOManager.loadClustersFromFile(inputFileName);
        if (!fullClustering) {
            logger.logAndStderr("Generating input statistics...");
            if (labels == null) {
                labels = getSortedLabels(FileIOManager.getAllSequences(inClusters));
            }
            if (!inGalaxy) {
                FileIOManager.saveInputStatistics(FileIOManager.getAllSequences(inClusters), labels, inputStatistics);
            }
        }
        List<Cluster> toCluster = new ArrayList<>();
        List<Cluster> other = new ArrayList<>();
        if (maxAlnLength == null) {
            maxAlnLength = setMaxAlnLength(inClusters);
            logger.logAndStderr("Maximal alignment length not set. Setting automatically to: " + maxAlnLength);
        }
        if (countThreshold == null) {
            if (sizeThreshold == null && partThreshold == null) {
                countThreshold = getDefaultCountThreshold(inClusters.size());
                logger.logAndStderr("Number of cluster cores not set. Setting automatically to: " + countThreshold);
            }
            if (sizeThreshold != null) {
                int count = 0;
                for (Cluster cl : inClusters) {
                    if (cl.size() >= sizeThreshold) {
                        count++;
                    }
                }
                countThreshold = count;
            }
            if (partThreshold != null) {
                countThreshold = (int) (partThreshold * inClusters.size());
            }
        }
        if (countThreshold >= inClusters.size()) {
            logger.logAndStderr("Error. 100% of initial clusters selected as cluster cores. This is not feasable. In most cases, "
                    + "this value should be less than 50%.");
            return;
        }
        if ((countThreshold * 2) >= inClusters.size()) {
            logger.logAndStderr("Warning. More than 50% of initial clusters selected as cluster cores. This may have negative influence or results");
        }
        Collections.sort(inClusters, Collections.reverseOrder(new ClusterSizeIdComparator()));
        toCluster.addAll(inClusters.subList(0, countThreshold));
        other.addAll(inClusters.subList(countThreshold, inClusters.size()));
        List<UniqueSequence> databaseSequences = new ArrayList<>();
        for (Cluster cl : other) {
            databaseSequences.addAll(cl.getSequences());
        }
        if (assignThresholdSequence == null) {
            assignThresholdSequence = setAssignThresholdSequence(inClusters);
            logger.logAndStderr("Assign threshold sequence not set. Setting automatically based on mean sequence length to: ");
            StringBuilder seq = new StringBuilder();
            for (double th : assignThresholdSequence) {
                seq.append(th).append(",");
            }
            logger.logAndStderr(seq.toString());
        }
        if (overlapThresholdSequence == null) {
            overlapThresholdSequence = setOverlapThresholdSequence(inClusters);
            StringBuilder seq = new StringBuilder();
            for (double th : overlapThresholdSequence) {
                seq.append(th).append(",");
            }
            logger.logAndStderr(seq.toString());
        }
        if (hhMergeThresholdSequence == null) {
            hhMergeThresholdSequence = setHhMergeThresholdSequence(inClusters);
            StringBuilder seq = new StringBuilder();
            for (double th : hhMergeThresholdSequence) {
                seq.append(th).append(",");
            }
            logger.logAndStderr(seq.toString());
        }

        if ((overlapThresholdSequence.length != assignThresholdSequence.length)
                || (hhMergeThresholdSequence.length != assignThresholdSequence.length)) {
            logger.logAndStderr("Error. Merge threshold sequence or overlap threshold "
                    + "sequence do not have the same length as assign threshold sequence.");
            return;
        }

        fullHHClustering = new boolean[overlapThresholdSequence.length];
        for (int i = 0; i < fullHHClustering.length; i++) {
            if (overlapThresholdSequence[i] == 0.0) {
                fullHHClustering[i] = true;
            } else {
                fullHHClustering[i] = false;
            }
        }

        int toBuild = 0;
        for (Cluster cl : toCluster) {
            if (!(cl.hasMSA())) {
                toBuild++;
            }
        }
        if (toBuild > 0) {
            logger.logAndStderr("Some cluster cores don't have MSAs. \nBuilding " + toBuild + " MSAs with Clustal...");
            ClustalRunner.mulitpleAlignment(toCluster);
        }
        Set<Integer> toReject = Hammock.getClusterIdsNotSatisfyingIc(toCluster);
        if (toReject.size() > 0) {
            logger.logAndStderr(toReject.size() + " clusters rejected because of match states and information content constraints.");
            List<Cluster> newList = new ArrayList<>();
            for (Cluster cl : toCluster) {
                if (toReject.contains(cl.getId())) {
                    databaseSequences.addAll(cl.getSequences());
                } else {
                    newList.add(cl);
                }
            }
            toCluster = newList;
        }
        if (toBuild > 0) {
            logger.logAndStderr("Saving newly aligned clusters...");
            if (!inGalaxy) {
                FileIOManager.saveClustersToFile(inClusters, initialAlignedCl);
                FileIOManager.saveClusterSequencesToCsv(inClusters, initialAlignedSequenceCsv, labels);
                FileIOManager.SaveClustersToCsv(inClusters, initialAlignedClusterCsv, labels);

            }
        }
        if (!inGalaxy) {
            for (Cluster cl : toCluster) {
                FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln",
                        initialAlnFolder + separatorChar + cl.getId() + ".aln");
            }
        }
        logger.logAndStderr("\nClustering in " + assignThresholdSequence.length + " rounds...");
        long time = System.currentTimeMillis();
        AssignmentResult result = IterativeHmmClusterer.iterativeHmmClustering(toCluster, databaseSequences, assignThresholdSequence, overlapThresholdSequence, hhMergeThresholdSequence, fullHHClustering, minMatchStates, minIc, Hammock.maxAlnLength, nThreads);
        List<Cluster> resultingClusters = result.getClusters();
        logger.logAndStderr("\nReady. Clustering time : " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Resulting clusers: " + resultingClusters.size());
        logger.logAndStderr("Containing " + getClustersUniqueSize(resultingClusters) + " unique sequences and " + getClustersTotalSize(resultingClusters) + " total sequences.");
        logger.logAndStderr("Unique sequences not assigned: " + result.getDatabaseSequences().size() + ", total sequences not assigned: " + Hammock.getTotalSequenceCount(result.getDatabaseSequences()));
        logger.logAndStderr("Saving results to outupt files...");
        if (!inGalaxy) {
            FileIOManager.saveClustersToFile(resultingClusters, finalCl);
            FileIOManager.saveUniqueSequencesToFasta(result.getDatabaseSequences(), finalRemainingSequences);
        }
        FileIOManager.saveClusterSequencesToCsv(resultingClusters, finalSequenceCsv, labels);
        FileIOManager.SaveClustersToCsv(resultingClusters, finalClustersCsv, labels);
        if (!inGalaxy) {
            for (Cluster cl : resultingClusters) {
                FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + separatorChar + cl.getId() + ".aln",
                        finalAlnFolder + separatorChar + cl.getId() + ".aln");
            }
        }
        logger.logAndStderr("Results in: " + finalSequenceCsv);
        logger.logAndStderr("and: " + finalClustersCsv);
        if (!inGalaxy) {
            logger.logAndStderr("and: " + finalCl);
        }
    }

    /**
     * Parses arguments common for all modes
     *
     * @param args
     */
    private static void parseCommonArgs(String[] args) {
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("-i") || args[i].equals("--input")) {
                if (args.length > i + 1) {
                    inputFileName = args[i + 1];
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-d") || args[i].equals("--outputDirectory")) {
                if (args.length > i + 1) {
                    workingDirectory = args[i + 1];
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-t") || args[i].equals("--threads")) {
                if (args.length > i + 1) {
                    nThreads = Integer.parseInt(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-l") || args[i].equals("--labels")) {
                if (args.length > i + 1) {
                    labelString = args[i + 1];
                    i++;
                    continue;
                }
            }

            if (args[i].equals("--galaxy")) {
                inGalaxy = true;
                continue;
            }

            if (args[i].equals("--goc")) {
                if (args.length > i + 1) {
                    galaxyFinalClustersCsv = args[i + 1];
                    i++;
                    continue;
                }
            }

            if (args[i].equals("--gos")) {
                if (args.length > i + 1) {
                    galaxyFinalSequenceCsv = args[i + 1];
                    i++;
                    continue;
                }
            }
        }
    }

    /**
     * Parses arguments special for greedy mode
     *
     * @param args
     */
    private static void parseGreedyArgs(String[] args) {
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("-f") || args[i].equals("--file_format")) {
                if (args.length > i + 1) {
                    inputType = args[i + 1];
                    i++;
                    continue;
                }
            }
            if (args[i].equals("-m") || args[i].equals("--matrix")) {
                if (args.length > i + 1) {
                    matrixFile = args[i + 1];
                    i++;
                    continue;
                }
            }
            if (args[i].equals("-g") || args[i].equals("--greedy_threshold")) {
                if (args.length > i + 1) {
                    greedyThreshold = Integer.decode(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-x") || args[i].equals("--max_shift")) {
                if (args.length > i + 1) {
                    maxShift = Integer.decode(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-p") || args[i].equals("--gap_penalty")) {
                if (args.length > i + 1) {
                    shiftPenalty = Integer.decode(args[i + 1]);
                    i++;
                }
            }
        }
    }

    /**
     * Parses arguments special for cluster mode
     *
     * @param args
     */
    private static void parseClusteringArgs(String[] args) {
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("-s") || args[i].equals("--size_threshold")) {
                if (args.length > i + 1) {
                    sizeThreshold = Integer.decode(args[i + 1]);
                    i++;
                    continue;
                }
            }
            if (args[i].equals("-c") || args[i].equals("--count_threshold")) {
                if (args.length > i + 1) {
                    countThreshold = Integer.decode(args[i + 1]);
                    i++;
                    continue;
                }
            }
            if (args[i].equals("-a") || args[i].equals("--part_threshold")) {
                if (args.length > i + 1) {
                    partThreshold = Double.parseDouble(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-n") || args[i].equals("--assign_thresholds")) {
                if (args.length > i + 1) {
                    assignThresholdSequenceString = args[i + 1].trim();
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-v") || args[i].equals("--overlap_thresholds")) {
                if (args.length > i + 1) {
                    overlapThresholdSequenceString = args[i + 1].trim();
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-r") || args[i].equals("--merge_thresholds")) {
                if (args.length > i + 1) {
                    hhMergeThresholdSequenceString = args[i + 1].trim();
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-b") || args[i].equals("--absolute_thresholds")) {
                relativeHHScore = false;
                relativeHmmScore = false;
            }

            if (args[i].equals("-e") || args[i].equals("--relative_thresholds")) {
                relativeHHScore = true;
                relativeHmmScore = true;
            }

            if (args[i].equals("-h") || args[i].equals("--min_match_states")) {
                if (args.length > i + 1) {
                    minMatchStates = Integer.decode(args[i + 1].trim());
                    if (minMatchStates <= 0) {
                        minMatchStates = 0;
                    }
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-y") || args[i].equals("--max_gap_proportion")) {
                if (args.length > i + 1) {
                    maxGapProportion = Double.parseDouble(args[i + 1].trim());
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-k") || args[i].equals("--min_ic")) {
                if (args.length > i + 1) {
                    minIc = Double.parseDouble(args[i + 1].trim());
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-j") || args[i].equals("--max_aln_length")) {
                if (args.length > i + 1) {
                    maxAlnLength = Integer.parseInt(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-u") || args[i].equals("--max_inner_gaps")) {
                if (args.length > i + 1) {
                    maxInnerGaps = Integer.parseInt(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-q") || args[i].equals("--extension_increase_length")) {
                extensionIncreaseLength = true;
            }

        }
    }

    /**
     * Checks if common arguments are sane and sets up file paths
     *
     */
    private static boolean checkCommonArgs() {
        if (inputFileName == null) {
            System.err.println("Error. Parameter input file (-i or --input) missing with no default.");
            return false;
        }

        if ((workingDirectory == null) && (!inGalaxy)) {
            int i = 1;
            String name = parentDir + separatorChar + "dist" + separatorChar + "Hammock_result_" + i;
            while (i < 9999) {
                name = parentDir + separatorChar + "dist" + separatorChar + "Hammock_result_" + i;
                File file = new File(name);
                if (!(file.exists())) {
                    file.mkdir();
                    break;
                }
                i++;
            }
            workingDirectory = name;
            System.err.println("Creating default output directory: " + name);
        }

        if (!inGalaxy) {
            logger = new Logger(workingDirectory + separatorChar + "run.log", false);
        } else {
            logger = new Logger(null, true);
        }

        if (labelString != null) {
            labels = new ArrayList<>();
            String[] splitLine = labelString.split(",");
            labels.addAll(Arrays.asList(splitLine));
        }

        initialAlignedCl = workingDirectory + separatorChar + "initial_clusters_aligned.cl";
        initialAlignedSequenceCsv = workingDirectory + separatorChar + "initial_clusters_aligned_sequences.csv";
        initialAlignedClusterCsv = workingDirectory + separatorChar + "initial_clusters.csv";
        finalCl = workingDirectory + separatorChar + "final_clusters.cl";
        finalRemainingSequences = workingDirectory + separatorChar + "final_remaining_sequences.fa";
        finalSequenceCsv = workingDirectory + separatorChar + "final_clusters_sequences.csv";
        finalClustersCsv = workingDirectory + separatorChar + "final_clusters.csv";
        inputStatistics = workingDirectory + separatorChar + "input_statistics.csv";

        threadPool = Executors.newFixedThreadPool(nThreads);
        return true;
    }

    /**
     * Checks if command line arguments for greedy mode are sane
     *
     */
    private static boolean checkGreedyArgs() {
        if (!((inputType.equals("fasta")) || (inputType.equals("seq")) || (inputType.equals("tab")))) {
            System.err.println("Error. Parameter -f value may be either \"fasta\", \"seq\" or \"tab\". No other values are allowed");
            return false;
        }
        return true;

    }

    /**
     * Checks if command line arguments for cluster mode are sane
     *
     */
    private static boolean checkClusteringArgs() {
        if (!checkEnvVariable()) {
            System.err.println("Error. Environmental variable \"HHLIB\" not set. Set it so that it contains path to: hhsuite_folder/lib/hh ");
            return false;
        }
        if ((partThreshold != null) && ((partThreshold > 1.0) || (partThreshold < 0.0))) {
            System.err.println("Error. Parameter -a (--part_threshold) must be within interval (0.0, 1.0)).");
            return false;
        }
        if (assignThresholdSequenceString != null) {
            assignThresholdSequence = parseThresholdSequence(assignThresholdSequenceString);
        }
        if (overlapThresholdSequenceString != null) {
            overlapThresholdSequence = parseThresholdSequence(overlapThresholdSequenceString);
        }
        if (hhMergeThresholdSequenceString != null) {
            hhMergeThresholdSequence = parseThresholdSequence(hhMergeThresholdSequenceString);
        }

        if (minIc > 4.3219280949) {
            System.err.println("Error. Minimal informationc content (-k) can't be more than 4.3219280949");
            return false;
        }

        if (!inGalaxy) {
            //Create directories for alignment files
            String name = workingDirectory + separatorChar + "alignments_initial";
            File file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
            initialAlnFolder = name;
            name = workingDirectory + separatorChar + "alignments_final";
            file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
            finalAlnFolder = name;

            name = workingDirectory + separatorChar + "alignments_other";
            file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
        }

        //delete contents of folders in temp
        FileIOManager.deleteFolderContents(Settings.getInstance().getFastaDirectory());
        FileIOManager.deleteFolderContents(Settings.getInstance().getHhDirectory());
        FileIOManager.deleteFolderContents(Settings.getInstance().getHhsearchOutDirectory());
        FileIOManager.deleteFolderContents(Settings.getInstance().getHmmDirectory());
        FileIOManager.deleteFolderContents(Settings.getInstance().getHmmsearchOutDirectory());
        FileIOManager.deleteFolderContents(Settings.getInstance().getMsaDirectory());

        return true;
    }

    /**
     * Sets default value for greedy clustering threshold
     *
     * @param sequences
     * @return
     */
    private static int setGreedyThreshold(Collection<UniqueSequence> sequences) {
        double meanLength = getMeanSequenceLength(sequences);
        int result = (int) Math.round(meanLength * 2.0);
        return result;
    }

    /**
     * Sets default value for maximal alignment length
     *
     * @param clusters
     * @return
     */
    private static int setMaxAlnLength(Collection<Cluster> clusters) {
        double meanLength = getClusterMeanSequenceLength(clusters);
        int result = (int) Math.round(meanLength * 2.0);
        return result;
    }

    /**
     * Sets default value for assign threshold sequence
     *
     * @param clusters
     * @return
     */
    private static double[] setAssignThresholdSequence(Collection<Cluster> clusters) {
        double meanLength = getClusterMeanSequenceLength(clusters);
        double[] result;
        if (relativeHmmScore) {
            result = new double[]{meanLength * 0.13, meanLength * 0.113, meanLength * 0.108};
        } else {
            result = new double[]{meanLength * 0.95, meanLength * 0.75, meanLength * 0.55};
        }
        for (int i = 0; i < result.length; i++) {
            result[i] = (double) Math.round(result[i] * 100) / 100;
        }
        return result;
    }

    /**
     * Sets default values for overlap threshold sequence
     *
     * @param clusters
     * @return
     */
    private static double[] setOverlapThresholdSequence(Collection<Cluster> clusters) {
        double[] result;
        if (assignThresholdSequence.length == 3) {
            logger.logAndStderr("Overlap threshold not set. Setting automatically based on mean sequence length to: ");
            double meanLength = getClusterMeanSequenceLength(clusters);
            if (relativeHHScore) {
                result = new double[]{meanLength * 0.09, meanLength * 0.075, 0.0};
            } else {
                result = new double[]{meanLength * 0.7, meanLength * 0.4, 0.0};
            }
            for (int i = 0; i < result.length; i++) {
                result[i] = (double) Math.round(result[i] * 100) / 100;
            }
        } else {
            logger.logAndStderr("Overlap threshold not set. Setting automatically based on assign threshold sequence to: ");
            result = new double[assignThresholdSequence.length];
            for (int i = 0; i < assignThresholdSequence.length; i++) {
                result[i] = assignThresholdSequence[i] * 0.75;
            }
            result[assignThresholdSequence.length - 1] = 0.0; //last zero for full clustering
        }
        return result;
    }

    /**
     * Sets default values for merge threshold sequence
     *
     * @param clusters
     * @return
     */
    private static double[] setHhMergeThresholdSequence(Collection<Cluster> clusters) {
        double[] result;
        if (assignThresholdSequence.length == 3) {
            logger.logAndStderr("Merge threshold not set. Setting automatically based on mean sequence length to: ");
            double meanLength = getClusterMeanSequenceLength(clusters);
            if (relativeHHScore) {
                result = new double[]{meanLength * 0.125, meanLength * 0.115, meanLength * 0.110};
            } else {
                result = new double[]{meanLength * 0.95, meanLength * 0.75, meanLength * 0.55};
            }
            for (int i = 0; i < result.length; i++) {
                result[i] = (double) Math.round(result[i] * 100) / 100;
            }
        } else {
            logger.logAndStderr("Merge threshold not set. Setting automatically based assign threshold sequence to: ");
            result = new double[assignThresholdSequence.length];
            for (int i = 0; i < assignThresholdSequence.length; i++) {
                result[i] = assignThresholdSequence[i] * 1.0;
            }
        }
        return result;
    }

    /**
     * Calculates mean sequence length
     *
     * @param sequences
     * @return
     */
    private static double getMeanSequenceLength(Collection<UniqueSequence> sequences) {
        int count = 0;
        int lengthSum = 0;
        for (UniqueSequence seq : sequences) {
            count++;
            lengthSum += seq.getSequence().length;
        }
        double meanLength = ((double) lengthSum) / ((double) count);
        return meanLength;
    }

    /**
     * Calculates mean sequence length
     *
     * @param clusters
     * @return
     */
    private static double getClusterMeanSequenceLength(Collection<Cluster> clusters) {
        List<UniqueSequence> sequences = new ArrayList<>();
        for (Cluster cl : clusters) {
            sequences.addAll(cl.getSequences());
        }
        return getMeanSequenceLength(sequences);
    }

    /**
     * returns all sequence labels present in this collecton of sequences sorted
     * from the most represented to the least
     *
     * @param sequences Clusters from which to get sorted list of labels
     * @return
     */
    private static List<String> getSortedLabels(Collection<UniqueSequence> sequences) {
        Map<String, Integer> labelMap = new HashMap<>();
        for (UniqueSequence seq : sequences) {
            for (Map.Entry<String, Integer> entry : seq.getLabelsMap().entrySet()) {
                Integer count = labelMap.get(entry.getKey());
                if (count == null) {
                    count = 0;
                }
                count += entry.getValue();
                labelMap.put(entry.getKey(), count);
            }
        }
        Map<String, Integer> sortedMap = new TreeMap<>(new ValueComparator(labelMap));
        sortedMap.putAll(labelMap);
        List<String> result = new ArrayList<>();
        for (String label : sortedMap.keySet()) {
            result.add(label);
        }
        return result;
    }

    /**
     * Parses threshold sequence (assign, overlap or merge) from command line
     * arguments
     *
     * @param sequenceString
     * @return
     */
    private static double[] parseThresholdSequence(String sequenceString) {
        String[] splitLine = sequenceString.split(",");
        double[] result = new double[splitLine.length];
        for (int i = 0; i < splitLine.length; i++) {
            result[i] = Double.parseDouble(splitLine[i]);
        }
        return result;
    }

    /**
     * Computes default value for count threshold
     *
     * @param nClusters
     * @return
     */
    private static int getDefaultCountThreshold(int nClusters) {
        int nCores = (int) ((double) nClusters * 0.025);
        if (((double) nClusters * 0.025) > 250) {
            nCores = 250;
        }
        if (((double) nClusters * 0.025) < 25) {
            nCores = 25;
            if (((double) nClusters * 0.25) < 25) {
                nCores = (int) ((double) nClusters * 0.25);
            }
        }
        return nCores;
    }

    private static List<UniqueSequence> filterSequencesForLabels(Collection<UniqueSequence> sequences, List<String> labels) {
        List<UniqueSequence> result = new ArrayList<>();
        for (UniqueSequence seq : sequences) {
            Map<String, Integer> newLabelsMap = new HashMap<>();
            for (String label : labels) {
                if (seq.getLabelsMap().containsKey(label)) {
                    newLabelsMap.put(label, seq.getLabelsMap().get(label));
                }
            }
            if (!(newLabelsMap.isEmpty())) {
                result.add(new UniqueSequence(seq.getSequenceString(), newLabelsMap));
            }
        }
        return result;
    }

    private static boolean checkEnvVariable() {
        Map<String, String> env = System.getenv();
        return env.containsKey("HHLIB");
    }


    private static Set<Integer> getClusterIdsNotSatisfyingIc(Collection<Cluster> clusters) throws DataException, IOException {
        Set<Integer> result = new HashSet<>();
        for (Cluster cl : clusters) {
            if (!cl.hasMSA()) {
                throw new DataException("Error, can't check match states for clusters without msas");
            }
            if (!FileIOManager.checkMatchStatesAndIc(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln", minMatchStates, minIc, maxGapProportion)) {
                result.add(cl.getId());
            }
        }
        return result;
    }

    private static int getTotalSequenceCount(Collection<UniqueSequence> sequences) {
        int count = 0;
        for (UniqueSequence seq : sequences) {
            for (int val : seq.getLabelsMap().values()) {
                count += val;
            }
        }
        return count;
    }

    private static int getClustersTotalSize(Collection<Cluster> clusters) {
        int size = 0;
        for (Cluster cl : clusters) {
            size += cl.size();
        }
        return size;
    }

    private static int getClustersUniqueSize(Collection<Cluster> clusters) {
        int size = 0;
        for (Cluster cl : clusters) {
            size += cl.getUniqueSize();
        }
        return size;
    }
}

/**
 * Compares clusters based on size. If sizes are equal, compares based on id.
 *
 */
class ClusterSizeIdComparator implements Comparator<Cluster> {

    @Override
    public int compare(Cluster o1, Cluster o2) {
        if (o1.size() == o2.size()) {
            return o2.getId() - o1.getId();
        } else {
            return o1.size() - o2.size();
        }
    }
}

/**
 * Compares sequences according to sum of counts of all labels. In case 
 * of equality, compares alphabetically.
 * @author Adam Krejci
 */
class UniqueSequenceSizeAlphabeticComparator implements Comparator<UniqueSequence>{

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        int result = new SizeComparator().compare(o1, o2);
        if (result == 0){
            result = new UniqueSequenceAlphabeticComparator().compare(o1, o2);
        }
        return result;
    }  
}

/**
 * Compares sequences alphabetically
 * @author Adam Krejci
 */
class UniqueSequenceAlphabeticComparator implements Comparator<UniqueSequence>{

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        return o1.getSequenceString().compareTo(o2.getSequenceString());
    }   
}
