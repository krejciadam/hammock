/**
 * Main class. Parses command line arguments and runs appropriate routines
 */
package cz.krejciadam.hammock;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
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
import java.util.Random;
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

    public static final char SEPARATOR_CHAR = File.separatorChar; //system separator
    public static final String CSV_SEPARATOR = "\t";
    private static final String PARENT_DIR = (new File(Hammock.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getParentFile().getPath());

    //common
    private static final String VERSION = "1.1.4";
    private static List<UniqueSequence> initialSequences = null;
    private static String inputFileName = null;
    public static String workingDirectory = null;
    public static int nThreads = 4;
    private static List<String> labels = null;
    private static String matrixFile = PARENT_DIR + SEPARATOR_CHAR + "matrices" + SEPARATOR_CHAR + "blosum62.txt";
    public static int[][] scoringMatrix = null;

    private static String labelString = null;
    private static Boolean fullClustering = null;
    public static ExecutorService threadPool = null;
    public static boolean inGalaxy = false;
    public static Random random = null;

    //common - preset
    private static String initialClustersSequencesCsv;
    private static String initialClusters;
    private static String initialClustersSequencesOrderedCsv;
    private static String finalRemainingSequences;
    private static String finalSequenceCsv;
    private static String finalClustersCsv;
    private static String finalSequenceOrderedCsv;
    private static String initialAlnFolder;
    private static String inputStatistics;
    private static String finalAlnFolder;
    public static Logger logger;
    public static Logger errorLogger = new Logger(null, false); //Only write to stderr until initialized
    public static int seed = 42;
    public static final String COUNT_MATRIX_FILE = PARENT_DIR + SEPARATOR_CHAR + "settings" + SEPARATOR_CHAR + "misc" + SEPARATOR_CHAR + "blosum62.freq_rownorm";
    public static String empiricalProbabsFile = null;
    public static String tempDirectory = "/tmp";
    
    //clustering
    private static String additionalSequencesPath = null;
    
    //full
    private static boolean useGreedy = false;
    private static boolean useClinkage = false;

    //greedy
    private static String inputType = "fasta";
    private static Integer sequenceClusteringThreshold = null;
    private static int shiftPenalty = 0;
    private static Integer maxShift = null;
    private static String order = "size";
    private static Integer initialClustersLimit = null;

    
    //clinkage
    private static int cacheSizeLimit = 1;
    
    //hmm
    private static boolean unique = false;
    private static Integer sizeThreshold = null;
    private static Integer countThreshold = null;
    private static Double partThreshold = null;
    private static String assignThresholdSequenceString = null;
    private static String overlapThresholdSequenceString = null;
    private static String hhMergeThresholdSequenceString = null;
    private static double[] assignThresholdSequence = null;
    private static double[] overlapThresholdSequence = null;
    private static double[] mergeThresholdSequence = null;
    private static boolean[] fullHHClustering = null;
    public static boolean relativeHHScore = false;
    public static boolean relativeHmmScore = false;
    public static Integer minConservedPositions = null;
    public static Double minIc = 1.2;
    public static Double maxGapProportion = 0.2;
    public static Integer maxAlnLength = null;
    public static int maxInnerGaps = 0;
    public static Boolean innerGapsAllowed = null;
    public static boolean extensionIncreaseLength = false;
    public static Map<String, String> hhSuiteEnv = null;
    //hmm - filtering
    public static boolean filterBeforeAssignment = false;
    public static int sequenceAddThreshold = 12;
    public static int gapOpenPenalty = -5;
    public static int gapExtendPenalty = -1;
    
    //correlation
    public static double minCorrelation = -1.0;
    public static int minClusterSize = 0;
    public static int minClusterUniqueSize = 0;
    
    //compare
    private static String databaseFileName = null;
    private static String secondClustersFileName = null;
    private static String compareResultsCsv = null;

    //galaxy
    private static String galaxyFinalClustersCsv;
    private static String galaxyFinalSequenceCsv;
    private static String galaxyFinalSequenceCsvOrdered;
    
    public static List<String> getLabels() {
        return labels;
    }

    public static void main(String[] args) throws IOException, HammockException, InterruptedException, ExecutionException, Exception {
        
        Alphabet sdm12 = AlphabetSdm12.getInstance();
        List<UniqueSequence> seqs = new ArrayList<>();
        seqs.add(new UniqueSequence("ADDGWGGGGG"));
        seqs.add(new UniqueSequence("ADDGWWW"));
        SequenceScorer scorer = new KmerScorer(3, sdm12);
        System.out.println(scorer.sequenceScore(seqs.get(0), seqs.get(1)));
        
        
        try {
            parseArgs(args);
        } catch (CLIException e){
            errorLogger.logAndStderr("Error in command line arguments: " + e.getMessage());
        } catch (FileFormatException e) {
            errorLogger.logAndStderr("Error. Probably wrong input file format? Run with --help for a brief description of command line parameters. Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            errorLogger.logAndStderr(errors.toString());
        } catch (NullPointerException e) {
            errorLogger.logAndStderr("Error. Maybe wrong input file format? Run with --help for a brief description of command line parameters. Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            errorLogger.logAndStderr(errors.toString());
        } catch (DataException e) {
            errorLogger.logAndStderr("Error. Maybe wrong input file format or wrong set of labels? Run with --help for a brief description of command line parameters. Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            errorLogger.logAndStderr(errors.toString());
        } catch (Exception e) {
            errorLogger.logAndStderr("Error. Run with --help for a brief description of command line parameters. Trace: \n");
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            errorLogger.logAndStderr(errors.toString());
        } finally{
            if (threadPool != null) {
                threadPool.shutdown();
            }
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
    private static void parseArgs(String[] args) throws IOException, FileFormatException, HammockException, InterruptedException, ExecutionException, CLIException, Exception {
        if (args.length < 1) {
            throw new CLIException("First argument must be a command. Run \"java -jar Hammock.jar --help\" for more info.");
        }

        if (args[0].equals("--help") || args[0].equals("-help") || args[0].equals("help") || args[0].equals("-h") || args[0].equals("-H")) {
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
                throw new CLIException("Error. Can't run in cluster mode in Galaxy");
            }
            parseClusteringArgs(args);
            if ((checkCommonArgs()) && (checkClusteringArgs())) {
                logger.logWithTime("Program started in mode \"cluster\".");
                logger.logWithoutTime("Command-line arguments: " + "\n" + argsString.toString() + "\n");
                logCommonParams(logger);
                runClustering();
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }

        if (args[0].equals("greedy")) {
            fullClustering = false;
            parseCommonArgs(args);
            if (inGalaxy) {
                throw new CLIException("Error. Can't run in greedy mode in Galaxy");
            }
            parseGreedyArgs(args);
            if (checkCommonArgs() && checkGreedyOrClinkageArgs()) {
                logger.logWithTime("Program started in mode \"greedy\".");
                logger.logWithoutTime("Command-line arguments: " + "\n" + argsString.toString() + "\n");
                logCommonParams(logger);
                logGreedyParams(logger);
                List<UniqueSequence> sequences = loadInputSequences(inputFileName, inputType, logger, labels);
                runGreedyClustering(sequences);
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }
        
        if (args[0].equals("clinkage")){
            fullClustering = false;
            parseCommonArgs(args);
            if (inGalaxy) {
                throw new CLIException("Error. Can't run in clinkage mode in Galaxy");
            }
            parseClinkageArgs(args);
            if (checkCommonArgs() && checkGreedyOrClinkageArgs()) {
                logger.logWithTime("Program started in mode \"clinkage\".");
                logger.logWithoutTime("Command-line arguments: " + "\n" + argsString.toString() + "\n");
                logCommonParams(logger);
                logClinkageParams(logger);
                List<UniqueSequence> sequences = loadInputSequences(inputFileName, inputType, logger, labels);
                runClinkageClustering(sequences);
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }

        if (args[0].equals("full")) {
            fullClustering = true;
            parseCommonArgs(args);
            parseClusteringArgs(args);
            parseClinkageArgs(args);
            parseGreedyArgs(args);
            
            if (checkCommonArgs() && checkGreedyOrClinkageArgs() && checkClusteringArgs()) {
                logger.logWithTime("Program started in mode \"full\".");
                logger.logWithoutTime("Command-line arguments: " + "\n" + argsString.toString() + "\n");
                logCommonParams(logger);
                if (useGreedy){
                    logGreedyParams(logger);
                } else{
                    logClinkageParams(logger);
                }
                runFull();
                logger.logWithTime("Program successfully ended.");
            }
            return;
        }
        
        if (args[0].equals("compare")){
            parseCommonArgs(args);
            parseCompareArgs(args);
            parseClusteringArgs(args);
            checkCommonArgs();
            checkClusteringArgs();
            runCompare();
            return;
        }
        throw new CLIException("The first argument must be a mode. Run \"java -jar Hammpck.jar --help\" for more info.");
    }

    /**
     * Displays help in command line
     */
    private static void displayHelp() {
        System.err.println("HAMMOCK - a Hidden Markov Model based peptide sequence clustering tool. For more details, see manual.pdf.");
        System.err.println("\nProject's home: https://github.com/krejciadam/hammock");
        System.err.println("Bug reports: https://github.com/krejciadam/hammock/issues");
        System.err.println("\nVersion: " + VERSION + "\n");
        System.err.println("\n------MANUAL------\n");
        System.err.println("Synopsis: java -jar Hammock.jar <mode> <param1> <param2> <param3>...");
        System.err.println("\nModes: \n");
        System.err.println("full\ngreedy\nclinkage\ncluster");
        System.err.println("\n------parameters common for all modes------\n");
        System.err.println("-i, --input <file>\n\tA path to an input file\n");
        System.err.println("-d, --output_directory <directory>\n\tA directory to store all output files in\n");
        System.err.println("-t, --threads <int>\n\tThe number of threads to use\n");
        System.err.println("-l, --labels <str,str,str...>\n\tA list of sequence labels to use\n");
        System.err.println("--temp <directory> \n\t A directory to store the temporal files in. Default: \\tmp");
        System.err.println("\n------parameters specific for greedy and clinkege modes------\n");
        System.err.println("-f, --file_format <[fasta,tab]>\n\tThe file format of input file specified by -i\n");
        System.err.println("-m, --matrix <file>\n\tA path to a substitution matrix file\n");
        System.err.println("-g, --alignment_threshold, (--greedy_threshold) <int>\n\tMinimal score needed for a sequence to join a cluster\n");
        System.err.println("-x, --max_shift <int>\n\tMaximal sequence-sequence shift. A nonnegative int\n");
        System.err.println("-p, --gap_penalty <int>\n\tThe penalty for each position of the sequence-sequence shift. A nonpositive int\n");
        System.err.println("\n------parameters specific for greedy mode------\n");
        System.err.println("-R, --order [size, alphabetic, random]\n\tThe order of sequences during greedy clustering\n");
        System.err.println("-S, --seed <int>\n\tA seed to make random processes deterministic (if -R random is in use)\n");
        System.err.println("--initial_clusters_limit <int>\n\tThe max. number of clusters resulting from gredy clustering\n");
        System.err.println("\n------parameters specific for clinkage mode------\n");
        System.err.println("-L, --cache_size_limit <int>\n\tIncrease for lower memory requirements (but lower speed as well) \n");
        System.err.println("\n------parameters specific for cluster mode------\n");
        System.err.println("-as, --additional_sequences\n\tA fasta file with sequences to add to the sequence pool\n");
        System.err.println("-U, --unique\n\tCluster cores will be selected on the basis of unique size instead of total size\n");
        //System.err.println("-a, --part_threshold <float [0, 1]>\n\tThe proportion of clusters to be used as cores\n");
        //System.err.println("-s, --size_threshold <int\n\tMinimal size (in unique sequences) of a cluster to be used a core\n");
        System.err.println("-c, --count_threshold <int>\n\tThis many largest (in terms of total or unique size) initial clusters will be used as cluster cores. \n");
        System.err.println("-n, --assign_thresholds <float,float,float...>\n\tA sequence of threshold values for sequence to cluster assignment\n");
        System.err.println("-v, --overlap thresholds <float,float,float...\n\tA sequence of threshold values of min. cluster overlap\n");
        System.err.println("-r, --merge_thresholds <float,float,float...>\n\tA sequence of threshold values for cluster-cluster comparisons\n");
        System.err.println("-e, --relative_thresholds\n\tAll thresholds are expressed as relative to the number of cluster match states\n");
        System.err.println("-b, --absolute_thresholds\n\tAll thresholds are expressed as absolute values\n");
        System.err.println("-y, --max_gap_proportion <float [0, 1]>\n\tMaximal proportion of gaps allowed in a match state/conserved MSA position\n");
        System.err.println("-k, --min_ic <float>\n\tMinimal information content allowed in a match state/conserved MSA positon\n");
        System.err.println("-h, --min_conserved_positions <int>\n\tClusters are not allowed to have less conserved MSA position. A conserved position is by -k, --min_ic and -y, --max_gap_proportion parameters.\n");
        System.err.println("-j, --max_aln_length <int>\n\tMaximal length of the cluster MSA including gaps\n");
        System.err.println("-u, --max inner_gaps <int>\n\tMaximal number of inner (non-trailing) gaps in the cluster MSA\n");
        System.err.println("-q, --extension_increase length>\n\tCluster extension step is allowed to increase MSA length\n");
        System.err.println("-C, --min_correlation <float [-1, 1]>   minimal correlation of cluster/cluster or cluster/sequence label count profiles\n");
        System.err.println("-M, --min_cluster_size <int>   clusters not reaching this size will be ommited from the result (and KLD calculation)\n");
        System.err.println("-N, --min_cluster_unique_size <int>   clusters not reaching this number of unique sequences will be ommited from the result (and KLD calculation)\n");
        System.err.println("\n------parameters specific for full mode------\n");
        System.err.println("--use_clinkage\n\tForce the use of clinkage clustering (for > 10 000 sequences)\n");
        System.err.println("--use_greedy\n\tForce the use of greedy clustering (for <= 10 000 sequences)\n");
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
            initialClustersSequencesCsv = Settings.getInstance().getTempDirectory() + SEPARATOR_CHAR + "greedy_aligned_cl";
            finalClustersCsv = galaxyFinalClustersCsv;
            finalSequenceCsv = galaxyFinalSequenceCsv;
            finalSequenceOrderedCsv = galaxyFinalSequenceCsvOrdered;

        }
        List<UniqueSequence> sequences = loadInputSequences(inputFileName, inputType, logger, labels);
        if (useGreedy){
            logger.logAndStderr("Greedy clustering forced.");
            runGreedyClustering(sequences);
        }
        else if (useClinkage){
            logger.logAndStderr("Clinkage clustering forced.");
            runClinkageClustering(sequences);
        } else if (sequences.size() <= 10000){
            logger.logAndStderr("Up to 10 000 unique sequences. Using clinkage clustering. Use --use_greedy to force greedy clustering");
            runClinkageClustering(sequences);
        } else{
            logger.logAndStderr("More than 10 000 unique sequences. Using greedy clustering. Use --use_clinkage to force clinkage clustering");
            runGreedyClustering(sequences);
        }
        inputFileName = initialClustersSequencesCsv;
        runClustering();
    }
    
    
    /**
     * Runs greedy clustering routine
     * @param sequences
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws DataException
     * @throws Exception 
     */
    private static void runGreedyClustering(List<UniqueSequence> sequences) throws IOException, InterruptedException, ExecutionException, DataException, Exception{
        prepareSequenceClustering(sequences);
        if (sequenceClusteringThreshold == null) {
            sequenceClusteringThreshold = setGreedyThreshold(sequences);
            logger.logAndStderr("Greedy clustering threshold not set. Setting automatically to: " + sequenceClusteringThreshold);
        }
        if (initialClustersLimit == null){
            initialClustersLimit = (int)Math.round(sequences.size() * 0.025);
            logger.logAndStderr("Initial greedy clusters limit not set. Setting automatically to: " + initialClustersLimit);
        }
        AligningSequenceScorer scorer = new ShiftedScorer(scoringMatrix, shiftPenalty, maxShift);
        SequenceClusterer clusterer = new LimitedGreedySequenceClusterer(scorer, sequenceClusteringThreshold, initialClustersLimit);
        
        logger.logAndStderr("Greedy clustering...");
        Long time = System.currentTimeMillis();
        sequences = UniqueSequence.sortSequences(sequences, order);
        
        List<Cluster> clusters = clusterer.cluster(sequences);
        
        logger.logAndStderr("Ready. Clustering time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Resulting clusers: " + clusters.size());
        logger.logAndStderr("Building MSAs...");
        List<Cluster> toBuildAlns = new ArrayList<>();
        List<Cluster> sizeOne = new ArrayList<>();
        for (Cluster cl : clusters){
            if (cl.getUniqueSize() > 1){
                toBuildAlns.add(cl);
            } else{
                sizeOne.add(cl);
            }
        }
        ClustalRunner.mulitpleAlignment(toBuildAlns);
        for (Cluster cl : sizeOne){
            cl.setAsHasMSA();
        }
        logger.logAndStderr("Ready. Total time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Saving results to output files...");
        FileIOManager.saveClusterSequencesToCsv(clusters, initialClustersSequencesCsv, labels);
        if (!inGalaxy) {
            FileIOManager.saveClusterSequencesToCsvOrdered(clusters, initialClustersSequencesOrderedCsv, labels, initialSequences);
            FileIOManager.SaveClustersToCsv(clusters, initialClusters, labels);
        }
        logger.logAndStderr("Greedy clustering results in: " + initialClusters);
        logger.logAndStderr("and: " + initialClustersSequencesCsv);
        logger.logAndStderr("and: " + initialClustersSequencesOrderedCsv);
    }
    
    
    /**
     * Runs complete linkage clustering routine
     * @param sequences
     * @throws HammockException
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws Exception 
     */
    private static void runClinkageClustering(List<UniqueSequence> sequences) throws HammockException, IOException, InterruptedException, ExecutionException, Exception{
        prepareSequenceClustering(sequences);
        
        if (sequenceClusteringThreshold == null) {
            sequenceClusteringThreshold = setClinkageThreshold(sequences);
            logger.logAndStderr("Clinkage clustering threshold not set. Setting automatically to: " + sequenceClusteringThreshold);
        }
        
        SequenceClusterer clusterer;
        ShiftedScorer scorer = new ShiftedScorer(scoringMatrix, shiftPenalty, maxShift);
        clusterer = new ClinkageSequenceClusterer(scorer, sequenceClusteringThreshold);
        logger.logAndStderr("Clinkage clustering...");
        Long time = System.currentTimeMillis();
        List<Cluster> clusters = clusterer.cluster(sequences);
        logger.logAndStderr("Ready. Clustering time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Resulting clusers: " + clusters.size());
        logger.logAndStderr("Building MSAs...");
        List<Cluster> toBuildAlns = new ArrayList<>();
        List<Cluster> sizeOne = new ArrayList<>();
        for (Cluster cl : clusters){
            if (cl.getUniqueSize() > 1){
                toBuildAlns.add(cl);
            } else{
                sizeOne.add(cl);
            }
        }
        ClustalRunner.mulitpleAlignment(toBuildAlns);
        for (Cluster cl : sizeOne){
            cl.setAsHasMSA();
        }
        logger.logAndStderr("Ready. Total time: " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Saving results to output files...");
        FileIOManager.saveClusterSequencesToCsv(clusters, initialClustersSequencesCsv, labels);
        if (!inGalaxy) {
            FileIOManager.saveClusterSequencesToCsvOrdered(clusters, initialClustersSequencesOrderedCsv, labels, initialSequences);
            FileIOManager.SaveClustersToCsv(clusters, initialClusters, labels);
        }
        logger.logAndStderr("Clinkage clustering results in: " + initialClusters);
        logger.logAndStderr("and: " + initialClustersSequencesCsv);
        logger.logAndStderr("and: " + initialClustersSequencesOrderedCsv);
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
        List<Cluster> inClusters = FileIOManager.loadClustersFromCsv(inputFileName, false);
        List<UniqueSequence> databaseSequences = new ArrayList<>();
        if (additionalSequencesPath != null){
            logger.logAndStderr("Loading additional sequences...");
            databaseSequences.addAll(FileIOManager.loadUniqueSequencesFromFasta(additionalSequencesPath));
        }
        if (!fullClustering) {
            logger.logAndStderr("Generating input statistics...");
            if (labels == null) {
                labels = getSortedLabels(FileIOManager.getAllSequences(inClusters));
            } else {
                inClusters = filterClustersForLabels(inClusters, labels);
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
        if (minConservedPositions == null) {
            minConservedPositions = setMinConservedPositions(inClusters);
            logger.logAndStderr("Minimal number of match states not set. Setting automatically to: " + minConservedPositions);
        }
        if (countThreshold == null) {
            if (sizeThreshold == null && partThreshold == null) {
                countThreshold = getDefaultCountThreshold(inClusters.size());
                logger.logAndStderr("Cluster cores count not set. Setting automatically to: " + countThreshold);
            }
            if (sizeThreshold != null) {
                int count = 0;
                for (Cluster cl : inClusters) {
                    if (unique){
                        if (cl.getUniqueSize() >= sizeThreshold){
                            count ++;
                        }
                    } else{
                        if (cl.size() >= sizeThreshold) {
                            count++;
                        }
                    }
                }
                countThreshold = count;
            }
            if (partThreshold != null) {
                countThreshold = (int) (partThreshold * inClusters.size());
            }
        }
        if (unique){
            Collections.sort(inClusters, Collections.reverseOrder(new ClusterUniqueSizeIdComparator()));
        } else{
            Collections.sort(inClusters, Collections.reverseOrder(new ClusterSizeIdComparator()));
        }
        toCluster.addAll(inClusters.subList(0, countThreshold));
        toCluster = FileIOManager.loadClusterAlignmentsFromFile(toCluster, inputFileName);
        other.addAll(inClusters.subList(countThreshold, inClusters.size()));
        for (Cluster cl : other) {
            databaseSequences.addAll(cl.getSequences());
        }
        if (assignThresholdSequence == null) {
            assignThresholdSequence = setAssignThresholdSequence(inClusters);
        }
        if (overlapThresholdSequence == null) {
            overlapThresholdSequence = setOverlapThresholdSequence(inClusters);
        }
        if (mergeThresholdSequence == null) {
            mergeThresholdSequence = setHhMergeThresholdSequence(inClusters);
        }

        if ((overlapThresholdSequence.length != assignThresholdSequence.length)
                || (mergeThresholdSequence.length != assignThresholdSequence.length)) {
            throw new CLIException("Error. Merge threshold sequence or overlap threshold "
                    + "sequence do not have the same length as assign threshold sequence.");
        }

        fullHHClustering = new boolean[overlapThresholdSequence.length];
        for (int i = 0; i < fullHHClustering.length; i++) {
            fullHHClustering[i] = overlapThresholdSequence[i] == 0.0;
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
                FileIOManager.saveClusterSequencesToCsv(inClusters, initialClustersSequencesCsv, labels);
                FileIOManager.SaveClustersToCsv(inClusters, initialClusters, labels);

            }
        }
        if (!inGalaxy) {
            for (Cluster cl : toCluster) {
                FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln",
                        initialAlnFolder + SEPARATOR_CHAR + cl.getId() + ".aln");
            }
        }

        logClusteringParams(logger);

        logger.logAndStderr("\nClustering in " + assignThresholdSequence.length + " rounds...");
        long time = System.currentTimeMillis();
        SequenceScorer scorer = null;
        if (Hammock.filterBeforeAssignment) {
            scorer = new LocalAlignmentScorer(scoringMatrix, gapOpenPenalty, gapExtendPenalty);
        }
        AssignmentResult result = IterativeHmmClusterer.iterativeHmmClustering(toCluster, databaseSequences, assignThresholdSequence, overlapThresholdSequence, mergeThresholdSequence, fullHHClustering, scorer, nThreads);
        List<Cluster> resultingClusters = result.getClusters();
        int origSize = resultingClusters.size();
        resultingClusters = filterClustersForSize(resultingClusters, minClusterSize, minClusterUniqueSize);
        if ((origSize - resultingClusters.size()) > 0){
            logger.logAndStderr("\n" + (origSize - resultingClusters.size()) + " clusters filtered out because of minimal size limits.");
        }
        logger.logAndStderr("\nReady. Clustering time : " + (System.currentTimeMillis() - time));
        logger.logAndStderr("Resulting clusers: " + resultingClusters.size());
        logger.logAndStderr("Containing " + getClustersUniqueSize(resultingClusters) + " unique sequences and " + getClustersTotalSize(resultingClusters) + " total sequences.");
        logger.logAndStderr("Unique sequences not assigned: " + result.getDatabaseSequences().size() + ", total sequences not assigned: " + Hammock.getTotalSequenceCount(result.getDatabaseSequences()));
        logger.logAndStderr("Saving results to outupt files...");
        if (!inGalaxy) {
            FileIOManager.saveUniqueSequencesToFasta(result.getDatabaseSequences(), finalRemainingSequences);
        }
        FileIOManager.saveClusterSequencesToCsv(resultingClusters, finalSequenceCsv, labels);
        FileIOManager.SaveClustersToCsv(resultingClusters, finalClustersCsv, labels);
        if (initialSequences != null) {
            FileIOManager.saveClusterSequencesToCsvOrdered(resultingClusters, finalSequenceOrderedCsv, labels, initialSequences);
        }
        if (!inGalaxy) {
            for (Cluster cl : resultingClusters) {
                FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + SEPARATOR_CHAR + cl.getId() + ".aln",
                        finalAlnFolder + SEPARATOR_CHAR + cl.getId() + ".aln");
            }
        }
        logger.logAndStderr("Results in: " + finalSequenceCsv);
        logger.logAndStderr("and: " + finalClustersCsv);
        if (initialSequences != null) {
            logger.logAndStderr("and: " + finalSequenceOrderedCsv);
        }
        if (!inGalaxy) {
            logger.logAndStderr("\nCalculating KLD...");
            List<String> paths = new ArrayList<>();
            int omitted = 0;
            for (Cluster cl : resultingClusters){
                if (cl.getUniqueSize() > 1){
                    paths.add(finalAlnFolder + SEPARATOR_CHAR + cl.getId() + ".aln");
                } else{
                    omitted += 1;
                }
            }
            if (omitted > 0){
                logger.logAndStderr(omitted + " clusters omitted from KLD calculation because each of them only contains a single unique sequence.");
            }
            double kld1 = Statistics.getMeanSystemKld(paths, false);
            double kld2 = Statistics.getMeanSystemKld(paths, true);
            logger.logAndStderr("Final system KLD over match state MSA positions: " + kld1);
            logger.logAndStderr("Final system KLD over all MSA positions: " + kld2);
        }
    }
    
    /**
     * Runs the compare routine
     * @throws IOException
     * @throws FileFormatException
     * @throws ExecutionException
     * @throws Exception 
     */
    private static void runCompare() throws IOException, FileFormatException, ExecutionException, Exception {
        if (secondClustersFileName != null){
            logger.logAndStderr("\nLoading cluster set 1...");
            List<Cluster> inClusters1 = FileIOManager.loadClustersFromCsv(inputFileName, true);
            ensureClusterMsas(inClusters1);
            logger.logAndStderr("\nLoading cluster set 2...");
            List<Cluster> inClusters2 = FileIOManager.loadClustersFromCsv(secondClustersFileName, true);
            ensureClusterMsas(inClusters2);
            logger.logAndStderr("\nBuilding HHs....");
            HHsuiteRunner.buildHHs(inClusters1, threadPool);
            HHsuiteRunner.buildHHs(inClusters2, threadPool);
            logger.logAndStderr("\nCalculating scores....");
            List<HHalignHit> hits = HHsuiteRunner.alignAllVsAll(inClusters1, inClusters2, threadPool);
            logger.logAndStderr("\nWriting results...");
            FileIOManager.saveHHAlignHitsToCsv(hits, inClusters1, inClusters2, compareResultsCsv);
        }
        
        else{
            logger.logAndStderr("\nLoading clusters...");
            List<Cluster> inClusters = FileIOManager.loadClustersFromCsv(inputFileName, true);
            logger.logAndStderr("Loading input sequences...");
            List<UniqueSequence> sequences = loadInputSequences(databaseFileName, inputType, logger, labels);
            //List<UniqueSequence> sequences = FileIOManager.loadUniqueSequencesFromFasta(databaseFileName);
        
            ensureClusterMsas(inClusters);
            HmmerRunner.buildHmms(inClusters);
            List<HmmsearchSequenceHit> hits = HmmerRunner.searchWithHmms(inClusters, sequences);
            FileIOManager.saveHmmsearchHitsToCsv(hits, compareResultsCsv, empiricalProbabsFile, inClusters.size(), sequences.size());
        }
        logger.logAndStderr("\nReady.");
        
    }
    
    /**
     * Loads unique sequences from an input file
     * @param inputFileName
     * @return
     * @throws FileFormatException
     * @throws HammockException
     * @throws IOException 
     */
    private static List<UniqueSequence> loadInputSequences(String inputFileName, String inputType, Logger logger, List<String> labels) throws FileFormatException, HammockException, IOException{
        List<UniqueSequence> sequences = null;
        logger.logAndStderr("Loading input sequences...");
        if (inputType.equals("fasta")) {
            sequences = FileIOManager.loadUniqueSequencesFromFasta(inputFileName);
        }
        if (inputType.equals("tab")) {
            sequences = FileIOManager.loadUniqueSequencesFromTable(inputFileName);
        }

        if (sequences == null) {
            throw new HammockException("Error, this should have been checked.");
        }
        logger.logAndStderr(sequences.size() + " unique sequences loaded.");
        logger.logAndStderr(Hammock.getTotalSequenceCount(sequences) + " total sequences loaded.");

        if (labels != null) {
            sequences = filterSequencesForLabels(sequences, labels);
        }

        logger.logAndStderr(sequences.size() + " unique sequences after non-specified labels filtered out");
        logger.logAndStderr(new Cluster(sequences, -1).size() + " total sequences after non-specified labels fileterd out");
        int minLength = Integer.MAX_VALUE; 
        int maxLength = Integer.MIN_VALUE;
        for (UniqueSequence seq : sequences){
            if (seq.getSequence().length > maxLength){
                maxLength = seq.getSequence().length;
            }
            if (seq.getSequence().length < minLength){
                minLength = seq.getSequence().length;
            }
        }
        logger.logAndStderr("Shortest sequence: " + minLength + " AA. Longest sequence: " + maxLength + " AA.");

        if (sequences.isEmpty()) {
            throw new FileFormatException("Error. No sequences (with specified labels) to cluster.");
        }
        return(sequences);
    }
    
    
    /**
     * Performs checks and prepares global variables for sequence (= greedy or clinkage) clustering
     * @param sequences
     * @throws IOException 
     */
    private static void prepareSequenceClustering(List<UniqueSequence> sequences) throws IOException{
        if (labels == null) {
            labels = getSortedLabels(sequences);
        }
        
        initialSequences = new ArrayList<>();
        initialSequences.addAll(sequences);

        if (maxShift == null) {
            maxShift = getMaxShift(sequences);
            logger.logAndStderr("Max shift not set. Setting automatically to: " + maxShift);
        } else {
            int correctMaxShift = checkMaxShift(sequences, maxShift);
            if (maxShift != correctMaxShift) {
                maxShift = correctMaxShift;
                logger.logAndStderr("Setting max shift to " + correctMaxShift + " as the length of the shortest sequence is only " + (correctMaxShift + 1));
            }
        }
        logger.logAndStderr("Generating input statistics...");
        if (!inGalaxy) {
            FileIOManager.saveInputStatistics(sequences, labels, inputStatistics);
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

            if (args[i].equals("-m") || args[i].equals("--matrix")) {
                if (args.length > i + 1) {
                    matrixFile = args[i + 1];
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

            if (args[i].equals("--goo")) {
                if (args.length > i + 1) {
                    galaxyFinalSequenceCsvOrdered = args[i + 1];
                    i++;
                    continue;
                }
            }
            
            if (args[i].equals("--use_greedy")){
                useGreedy = true;
                continue;
            }
            
            if (args[i].equals("--use_clinkage")){
                useClinkage = true;
            }  
            
            if (args[i].equals("--temp")){
                tempDirectory = args[i + 1];
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
            if (args[i].equals("-g") || args[i].equals("--greedy_threshold") || args[i].equals("--alignment_threshold")) {
                if (args.length > i + 1) {
                    sequenceClusteringThreshold = Integer.decode(args[i + 1]);
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

            if (args[i].equals("-R") || args[i].equals("--order")) {
                if (args.length > i + 1) {
                    order = args[i + 1];
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-S") || args[i].equals("--seed")) {
                if (args.length > i + 1) {
                    seed = Integer.decode(args[i + 1]);
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
            
            if (args[i].equals("--initial_clusters_limit")) {
                if (args.length > i + 1) {
                    initialClustersLimit = Integer.decode(args[i + 1]);
                    i++;
                }
            }
        }
    }
    
    private static void parseClinkageArgs(String[] args){
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("-f") || args[i].equals("--file_format")) {
                if (args.length > i + 1) {
                    inputType = args[i + 1];
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
            
            if (args[i].equals("-g") || args[i].equals("--greedy_threshold") || args[i].equals("--alignment_threshold")) {
                if (args.length > i + 1) {
                    sequenceClusteringThreshold = Integer.decode(args[i + 1]);
                    i++;
                }
            }
            
            if (args[i].equals("-L") || args[i].equals("--cache_size_limit")) {
                if (args.length > i + 1) {
                    cacheSizeLimit = Integer.decode(args[i + 1]);
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
            
            if (args[i].equals("-as") || args[i].equals("--additional_sequences")) {
                if (args.length > i + 1) {
                    additionalSequencesPath = args[i + 1];
                    i++;
                    continue;
                }
            }
            
            if (args[i].equals("-U") || args[i].equals("--unique")) {
                unique = true;
                continue;
            }
            
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

            if (args[i].equals("-h") || args[i].equals("--min_conserved_positions")) {
                if (args.length > i + 1) {
                    minConservedPositions = Integer.decode(args[i + 1].trim());
                    if (minConservedPositions <= 0) {
                        minConservedPositions = 0;
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

            if (args[i].equals("-C") || args[i].equals("--min_correlation")) {
                if (args.length > i + 1) {
                    minCorrelation = Double.parseDouble(args[i + 1]);
                    i++;
                    continue;
                }
            }

            if (args[i].equals("-q") || args[i].equals("--extension_increase_length")) {
                extensionIncreaseLength = true;
                continue;
            }
            
            if (args[i].equals("-M") || args[i].equals("--min_cluster_size")) {
                if (args.length > i + 1) {
                    minClusterSize = Integer.parseInt(args[i + 1]);
                    i++;
                    continue;
                }
            }
            
            if (args[i].equals("-N") || args[i].equals("--min_cluster_unique_size")) {
                if (args.length > i + 1) {
                    minClusterUniqueSize = Integer.parseInt(args[i + 1]);
                    i++;
                }
            }

        }
    }
    
    private static void parseCompareArgs(String[] args) {
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("-f") || args[i].equals("--file_format")) {
                if (args.length > i + 1) {
                    inputType = args[i + 1];
                    i++;
                    continue;
                }
            } 
            if (args[i].equals("-D") || args[i].equals("--database_file")) {
                if (args.length > i + 1) {
                    databaseFileName = args[i + 1];
                    i++;
                    continue;
                }
            }
            if (args[i].equals("-E") || args[i].equals("--empirical_probabs_file")) {
                if (args.length > i + 1) {
                    empiricalProbabsFile = args[i + 1];
                    i++;
                }
            }
            if (args[i].equals("-i2") || args[i].equals("--input2")) {
                if (args.length > i + 1) {
                    secondClustersFileName = args[i + 1];
                    i++;
                }
            }
        }
    }
    /**
     * Checks if common arguments are sane and sets up file paths
     *
     */
    private static boolean checkCommonArgs() throws IOException, HammockException {
        if (inputFileName == null) {
            throw new CLIException("Error. Parameter input file (-i or --input) missing with no default.");
        }
        
        if ((workingDirectory != null) && (!inGalaxy)){
            File file = new File(workingDirectory);
            if (file.exists()){
                throw new CLIException("Error. Output directory exists. Exiting to prevent data loss.");
            }else{
                file.mkdir();
            }
        }

        if ((workingDirectory == null) && (!inGalaxy)) {
            int i = 1;
            String name = PARENT_DIR + SEPARATOR_CHAR + "dist" + SEPARATOR_CHAR + "Hammock_result_" + i;
            while (i < 9999) {
                name = PARENT_DIR + SEPARATOR_CHAR + "dist" + SEPARATOR_CHAR + "Hammock_result_" + i;
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
            logger = new Logger(workingDirectory + SEPARATOR_CHAR + "run.log", false);
            errorLogger = logger;
        } else {
            logger = new Logger(null, true); //errorLogger still writers to stderr
        }

        logger.logAndStderr("\nHammock version " + VERSION + " Run with --help for a brief description of command line parameters.\n");
        
        if (labelString != null) {
            labels = new ArrayList<>();
            String[] splitLine = labelString.split(",");
            labels.addAll(Arrays.asList(splitLine));
        }

        random = new Random(seed);
        initialClustersSequencesCsv = workingDirectory + SEPARATOR_CHAR + "initial_clusters_sequences.tsv";
        initialClustersSequencesOrderedCsv = workingDirectory + SEPARATOR_CHAR + "initial_clusters_sequences_original_order.tsv";
        initialClusters = workingDirectory + SEPARATOR_CHAR + "initial_clusters.tsv";
        finalRemainingSequences = workingDirectory + SEPARATOR_CHAR + "final_remaining_sequences.fa";
        finalSequenceCsv = workingDirectory + SEPARATOR_CHAR + "final_clusters_sequences.tsv";
        finalSequenceOrderedCsv = workingDirectory + SEPARATOR_CHAR + "final_clusters_sequences_original_order.tsv";
        finalClustersCsv = workingDirectory + SEPARATOR_CHAR + "final_clusters.tsv";
        compareResultsCsv = workingDirectory + SEPARATOR_CHAR + "comparison_results.tsv";
        inputStatistics = workingDirectory + SEPARATOR_CHAR + "input_statistics.tsv";

        threadPool = Executors.newFixedThreadPool(nThreads);
        scoringMatrix = FileIOManager.loadScoringMatrix(matrixFile);    
        return true;
    }

    /**
     * Checks if command line arguments for greedy mode are sane
     *
     */
    private static boolean checkGreedyOrClinkageArgs() throws CLIException {
        if (!((inputType.equals("fasta")) || (inputType.equals("seq")) || (inputType.equals("tab")))) {
            throw new CLIException("Error. Parameter -f value may be either \"fasta\", \"seq\" or \"tab\". No other values are allowed");
        }
        return true;
    }

    /**
     * Checks if command line arguments for cluster mode are sane
     *
     */
    private static boolean checkClusteringArgs() throws IOException, InterruptedException, CLIException {
        Map<String, String> env = System.getenv();
        if (!env.containsKey("HHLIB")){
            hhSuiteEnv = new HashMap<>();
            hhSuiteEnv.put("HHLIB", PARENT_DIR + SEPARATOR_CHAR + "hhsuite-2.0.16" + SEPARATOR_CHAR + "lib" + SEPARATOR_CHAR + "hh" + SEPARATOR_CHAR);
            logger.logAndStderr("HHLIB env. variable not set. Setting to: " + hhSuiteEnv.get("HHLIB"));
        }
        if (!inGalaxy) {
            checkExternalProgram(Settings.getInstance().getClustalCommand(), Arrays.asList("-h"), "Clustal Omega", null);
            checkExternalProgram(Settings.getInstance().getHmmbuildCommand(), Arrays.asList("-h"), "hmmbuild", null);
            checkExternalProgram(Settings.getInstance().getHmmsearchCommand(), Arrays.asList("-h"), "hmmsearch", null);
            checkExternalProgram(Settings.getInstance().getHhmakeCommand(), Arrays.asList("-h"), "hhmake", hhSuiteEnv);
            checkExternalProgram(Settings.getInstance().getHhsearchCommand(), Arrays.asList("-h"), "hhsearch", hhSuiteEnv);
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
            mergeThresholdSequence = parseThresholdSequence(hhMergeThresholdSequenceString);
        }

        if (minIc > 4.3219280949) {
            throw new CLIException("Error. Minimal informationc content (-k) can't be more than 4.3219280949");
        }

        if (!inGalaxy) {
            //Create directories for alignment files
            String name = workingDirectory + SEPARATOR_CHAR + "alignments_initial";
            File file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
            initialAlnFolder = name;
            name = workingDirectory + SEPARATOR_CHAR + "alignments_final";
            file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
            finalAlnFolder = name;

            name = workingDirectory + SEPARATOR_CHAR + "alignments_other";
            file = new File(name);
            if (!(file.exists())) {
                file.mkdir();
            }
        }
        
        if (maxInnerGaps > 0){
            innerGapsAllowed = true;
        } else{
            innerGapsAllowed = false;
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
    
    private static void ensureClusterMsas(List<Cluster> clusters) throws ExecutionException, Exception{
        int toBuild = 0;
        for (Cluster cl : clusters) {
            if (!(cl.hasMSA())) {
                toBuild++;
            }
        }
        if (toBuild > 0) {
            logger.logAndStderr("Some cluster cores don't have MSAs. \nBuilding " + toBuild + " MSAs with Clustal...");
            ClustalRunner.mulitpleAlignment(clusters);
        }
    }
    

    /**
     * Check if external program returns no errors. Returns true if the length
     * of stderr is zero, throws an exception otherwise.
     *
     * @param command A command to be run
     * @param args all command args
     * @param programName Program name to be displayed in the exception message
     * @return true, if stderr of the program is of length 0, throws an
     * exception otherwise.
     * @throws IOException if length of stderr is not 0.
     * @throws InterruptedException
     */
    private static boolean checkExternalProgram(String command, List<String> args, String programName,  Map<String, String> env) throws IOException, InterruptedException {
        ByteArrayOutputStream errorBaos = new ByteArrayOutputStream();
        PrintStream errorStream = new PrintStream(errorBaos);
        ByteArrayOutputStream outputBaos = new ByteArrayOutputStream();
        PrintStream outputStream = new PrintStream(outputBaos);

        try {
            ExternalProcessRunner.runProcess(command, args, outputStream, errorStream, "", env);
        } catch (IOException e) {
            throw new IOException("Error, can't run " + programName + ". Check, if the program is installed properly and and runnable.\n"
                    + "Current path to " + programName + ": " + Settings.getInstance().getClustalCommand() + "\n"
                    + "This path can be changed in settings.prop file.\n"
                    + "Original excepton message: \n"
                    + e.getMessage());
        }
        if (errorBaos.toString().length() == 0) {
            return (true);
        } else {
            throw new IOException("Error, " + programName + " won't run properly. Make sure it runs and try again. \n"
                    + "original error message: \n" + errorBaos.toString() + ""
                    + "original stdout message: \n" + outputBaos.toString());
        }
    }

    /**
     * Sets default value for greedy clustering threshold
     *
     * @param sequences
     * @return
     */
    private static int setGreedyThreshold(Collection<UniqueSequence> sequences) {
        double meanLength = getMeanSequenceLength(sequences);
        int result = (int) Math.round(meanLength * 1.7);
        return result;
    }
    
    private static int setClinkageThreshold(Collection<UniqueSequence> sequences){
        double meanLength = getMeanSequenceLength(sequences);
        int result = (int) Math.round(meanLength * 1.7);
        return result;
    }

    private static int checkMaxShift(Collection<UniqueSequence> sequences, int maxShift) {
        int minLength = Integer.MAX_VALUE;
        for (UniqueSequence seq : sequences) {
            minLength = Math.min(minLength, seq.getSequence().length);
        }
        return Math.min(maxShift, (minLength - 1));
    }

    private static int getMaxShift(Collection<UniqueSequence> sequences) {
        double meanLength = getMeanSequenceLength(sequences);
        int result = (int) Math.round(meanLength / 4);
        result = checkMaxShift(sequences, result);
        return (result);
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

    private static int setMinConservedPositions(Collection<Cluster> clusters) {
        double meanLength = getClusterMeanSequenceLength(clusters);
        int result = (int) Math.round(meanLength / 3);
        return result;
    }

    /**
     * Sets default value for assign threshold sequence
     *
     * @param clusters
     * @return
     */
    private static double[] setAssignThresholdSequence(Collection<Cluster> clusters) {
        logger.logAndStderr("Assign threshold not set. Setting automatically to: ");
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
        logger.logAndStderr(formatDoubleSequence(result).toString());
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
            logger.logAndStderr("Overlap threshold not set. Setting automatically to: ");
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
        logger.logAndStderr(formatDoubleSequence(result).toString());
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
            logger.logAndStderr("Merge threshold not set. Setting automatically based on average sequence length to: ");
            double meanLength = getClusterMeanSequenceLength(clusters);
            if (relativeHHScore) {
                result = new double[]{meanLength * 0.125, meanLength * 0.115, meanLength * 0.110};
            } else {
                result = new double[]{meanLength * 1.0, meanLength * 0.9, meanLength * 0.8};
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
        logger.logAndStderr(formatDoubleSequence(result).toString());
        return result;
    }
    
    private static StringBuilder formatDoubleSequence(double[] numbers){
        StringBuilder seq = new StringBuilder();
            for (int i = 0; i < numbers.length; i++) {
                seq.append(numbers[i]);
                if (i < numbers.length - 1){
                    seq.append(",");
                }
            }
          return(seq);
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

    private static List<Cluster> filterClustersForLabels(Collection<Cluster> clusters, List<String> labels) throws FileFormatException {
        List<Cluster> result = new ArrayList<>();
        for (Cluster cl : clusters) {
            result.add(new Cluster(filterSequencesForLabels(cl.getSequences(), labels), cl.getId()));
        }
        return (result);
    }
    
    private static List<Cluster> filterClustersForSize(Collection<Cluster> clusters, int minSize, int minUniqueSize){
        List<Cluster> result = new ArrayList<>();
        for (Cluster cl : clusters){
            if(cl.size() >= minSize && cl.getUniqueSize() >= minUniqueSize){
                result.add(cl);
            }
        }
        return(result);
    }

    private static List<UniqueSequence> filterSequencesForLabels(Collection<UniqueSequence> sequences, List<String> labels) throws FileFormatException {
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
    
    private static Set<Integer> getClusterIdsNotSatisfyingIc(Collection<Cluster> clusters) throws DataException, IOException {
        Set<Integer> result = new HashSet<>();
        for (Cluster cl : clusters) {
            if (!cl.hasMSA()) {
                throw new DataException("Error, can't check match states for clusters without msas");
            }
            if (!FileIOManager.checkConservedStates(FileIOManager.getAlignmentLines(cl), minConservedPositions, minIc, maxGapProportion)) {
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

    private static void logCommonParams(Logger logger) {
        StringBuilder message = new StringBuilder();
        message.append("\nComplete list of input/output parameters: \n");
        message.append("-i, --input ").append(inputFileName).append("\n");
        message.append("-d, --output_directory ").append(workingDirectory).append("\n");
        message.append("-t, --thread ").append(nThreads).append("\n");
        message.append("-l, --labels ").append(labels).append("\n");
        message.append("\n");
        logger.logWithoutTime(message.toString());
    }

    private static void logGreedyParams(Logger logger) {
        StringBuilder message = new StringBuilder();
        message.append("\nComplete list of greedy clustering parameters: \n");

        message.append("-f, --file_format ").append(inputType).append("\n");
        message.append("-m, --matrix ").append(matrixFile).append("\n");
        message.append("-g, --greedy_threshold ").append(sequenceClusteringThreshold).append("\n");
        message.append("-x, --max_shift ").append(maxShift).append("\n");
        message.append("-p, --gap_penalty ").append(shiftPenalty).append("\n");
        message.append("-R, --order ").append(order).append("\n");
        message.append("-S, --seed ").append(seed).append("\n");
        message.append("\n");
        logger.logWithoutTime(message.toString());
    }
    
    private static void logClinkageParams(Logger logger){
        StringBuilder message = new StringBuilder();
        message.append("\nComplete list of clinkage clustering parameters: \n");
        
        message.append("-f, --file_format ").append(inputType).append("\n");
        message.append("-m, --matrix ").append(matrixFile).append("\n");
        message.append("-g, --alignment_threshold (--greedy_threshold)").append(sequenceClusteringThreshold).append("\n");
        message.append("-x, --max_shift ").append(maxShift).append("\n");
        message.append("-p, --gap_penalty ").append(shiftPenalty).append("\n");
        message.append("-C, --cache_size_limit ").append(cacheSizeLimit).append("\n");
        
        message.append("\n");
        logger.logWithoutTime(message.toString());
    }

    private static void logClusteringParams(Logger logger) {
        StringBuilder message = new StringBuilder();
        message.append("\nComplete list of HMM-based clustering parameters: \n");
        message.append("-a, --part_threshold ").append(partThreshold).append("\n");
        message.append("-s, --size_threshold ").append(sizeThreshold).append("\n");
        message.append("-c, --count_threshold ").append(countThreshold).append("\n");
        message.append("-n, --assign_thresholds ");
        for (double th : assignThresholdSequence) {
            message.append(th).append(",");
        }
        message.append("\n");
        message.append("-v, --overlap_thresholds ");
        for (double th : overlapThresholdSequence) {
            message.append(th).append(",");
        }
        message.append("\n");
        message.append("-r, --merge_thresholds ");
        for (double th : mergeThresholdSequence) {
            message.append(th).append(",");
        }
        message.append("\n");
        message.append("-e, --relative_thresholds ").append(relativeHmmScore).append("\n");
        message.append("-b, --absolute_thresholds ").append(!relativeHmmScore).append("\n");
        message.append("-h, --min_conserved_positions ").append(minConservedPositions).append("\n");
        message.append("-y, --max_gap_proportion ").append(maxGapProportion).append("\n");
        message.append("-k, --min_ic ").append(minIc).append("\n");
        message.append("-j, --max_aln_length ").append(maxAlnLength).append("\n");
        message.append("-u, --max_inner_gaps ").append(maxInnerGaps).append("\n");
        message.append("-q, --extension_increase_length ").append(extensionIncreaseLength).append("\n");
        message.append("\n");
        logger.logWithoutTime(message.toString());
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
 * Compares clusters based on unique size. If sizes are equal, compares based on id.
 *
 */
class ClusterUniqueSizeIdComparator implements Comparator<Cluster> {

    @Override
    public int compare(Cluster o1, Cluster o2) {
        if (o1.getUniqueSize() == o2.getUniqueSize()) {
            return o2.getId() - o1.getId();
        } else {
            return o1.getUniqueSize() - o2.getUniqueSize();
        }
    }
}

