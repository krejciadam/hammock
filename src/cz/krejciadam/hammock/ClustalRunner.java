/*
 * Constructs multiple alignments running Clustal as external process. 
 */
package cz.krejciadam.hammock;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class ClustalRunner {

    //Refactor note: Should not modify clusters and return void. Shoud return a new collection of modified clusters
    /**
     * Generates multiple alignment files for a collection of clusters using
     * Clustal Omega. MSA files will be written in appropriate directory
     * specified in Settings. Sets hasMSA() to true to clusters when MSA
     * construction is successful. Input clusters must be modifiable.
     *
     * @param clusters Clusters for which MSAs will be generated
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static void mulitpleAlignment(Collection<Cluster> clusters) throws InterruptedException, ExecutionException, Exception {
        List<Callable<Void>> callables = new ArrayList<>();
        for (Cluster cl : clusters) {
            callables.add(new SingleThreadClustalRunner(cl));
        }
        ParallelRunner.ExecuteInParallel(callables, Hammock.threadPool);
    }

    /**
     * Generates multiple alignment file for a cluster using Clustal Omega. MSA
     * files will be written in appropriate directory specified in Settings.
     * Sets hasMSA after successfull completion. Cluster must be modifiable.
     *
     * @param cluster Cluster for which MSA will be generated
     * @throws IOException
     * @throws InterruptedException
     */
    public static void multipleAlignment(Cluster cluster) throws IOException, InterruptedException {
        SingleThreadClustalRunner runner = new SingleThreadClustalRunner(cluster);
        runner.call();
    }

    /**
     * into each cluster (key in extensionMap), inserts all appropriate
     * sequences (corresponding value in extension map). Runs clustal to generate
     * new MSA for every such extended cluster. Sequences are inserted into
     * clusters one by one, starting from the most similar ones. Sequences
     * causing a cluster to violate minMatchStates are rejected and returned as
     * rejected sequences in resulting object
     *
     * @param extensionMap 
     * @return
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static ExtendedClusters extendClusters(Map<Cluster, List<HmmsearchSequenceHit>> extensionMap, int minMatchStates, Double minIc, SequenceScorer scorer) throws InterruptedException, ExecutionException, DataException {
        List<Cluster> resultingClusters = new ArrayList<>();
        List<UniqueSequence> rejectedSequences = new ArrayList<>();
        CompletionService<ExtendClusterResult> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        for (Map.Entry<Cluster, List<HmmsearchSequenceHit>> entry : extensionMap.entrySet()) {
            List<UniqueSequence> sequences = new ArrayList<>();
            List<HmmsearchSequenceHit> hits = entry.getValue();
            Collections.sort(hits, Collections.reverseOrder());
            for (HmmsearchSequenceHit hit : hits) {
                sequences.add(hit.getSequence());
            }

            List<UniqueSequence> selectedSequences = new ArrayList<>();
            selectedSequences.add(sequences.get(0));
            if (sequences.size() >= 2) {
                if (Hammock.filterBeforeAssignment) {
                    for (int i = 1; i < sequences.size(); i++) {
                        if (scorer.sequenceScore(sequences.get(0), sequences.get(i)) >= Hammock.sequenceAddThreshold) {
                            selectedSequences.add(sequences.get(i));
                        } else {
                            rejectedSequences.add(sequences.get(i));
                        }
                    }

                } else {
                    selectedSequences.addAll(sequences.subList(1, sequences.size()));
                }
            }
            resultPool.submit(new SingleThreadExtendClusterRunnerClustal(entry.getKey(), selectedSequences, minMatchStates, minIc));
        }
        for (int i = 0; i < extensionMap.size(); i++) {
            ExtendClusterResult result = resultPool.take().get();
            resultingClusters.add(result.getCluster());
            rejectedSequences.addAll(result.getRejectedSequences());
        }
        return new ExtendedClusters(resultingClusters, rejectedSequences);
    }
}

/**
 * Runs single instance of Clustal as external process on single thread.
 *
 * @author Adam Krejci
 */
class SingleThreadClustalRunner implements Callable<Void> {

    private final Cluster cluster;

    public SingleThreadClustalRunner(Cluster cluster) {
        this.cluster = cluster;
    }

    /**
     * Method builds MSA for a single cluster, calls cluster's setAsHasMSA()
     * method if successfull. Cluster must be modifiable. File with resulting
     * MSA is saved in appropriate folder (according to Settings) with name in
     * format: [cluster ID].aln If the cluster already has MSA, nothing happens.
     * If the cluster only contains 1 uniqueSequence, its content is written as
     * fasta file in [cluster ID].aln file. If a cluster hasMSA() and a file
     * with appropriate name already exists, new alignment is constructed and
     * old file is overwritten without warning.
     *
     * @throws Exception
     */
    @Override
    public Void call() throws IOException, InterruptedException {
        if (cluster.hasMSA()) {  //no need to build MSA if we already have it
            return null;
        }
        if (cluster.getUniqueSize() <= 1) {
            FileIOManager.saveClusterToFasta(cluster,
                    Settings.getInstance().getMsaDirectory() + cluster.getId() + ".aln"); //1-sequence MSA is simply fasta with this sequence
            cluster.setAsHasMSA();
            return null;
        }
        String fastaPath = Settings.getInstance().getFastaDirectory() + cluster.getId() + ".fa";
        FileIOManager.saveClusterToFasta(cluster, fastaPath);
        List<String> parameters = new ArrayList<>();
        parameters.add("-i");
        parameters.add(fastaPath);
        parameters.add("-o");
        parameters.add(Settings.getInstance().getMsaDirectory() + cluster.getId() + ".aln");
        parameters.add("--force");
        List<String> clustalParameters = Settings.getInstance().getClustalParameters();
        if (clustalParameters != null) {
            parameters.addAll(clustalParameters);
        }
        ExternalProcessRunner.runProcess(Settings.getInstance().getClustalCommand(), parameters, System.out, System.err, "While builidng msa for cluster " + cluster.getId() + " ");
//        ExternalProcessRunner.runProcess(Settings.getInstance().getClustalCommand(), parameters, System.out, null, null);
        cluster.setAsHasMSA();
        return null;
    }
}


/**
 * Runs single instance of Clustal routine for cluster extension as external
 * process.
 *
 * @author akrejci
 */
class SingleThreadExtendClusterRunnerClustal implements Callable<ExtendClusterResult> {

    private final Cluster extendedCluster;
    private final List<UniqueSequence> insertedSequences;
    private final Integer minMatchStates;
    private final Double minIc;

    public SingleThreadExtendClusterRunnerClustal(Cluster extendedCluster, List<UniqueSequence> insertedSequences, Integer minMatchStates, Double minIc) {
        this.extendedCluster = extendedCluster;
        this.insertedSequences = insertedSequences;
        this.minMatchStates = minMatchStates;
        this.minIc = minIc;

    }

    /**
     * Method inserts all sequences provided into cluster provided and builds
     * MSA for this clsuter. The process is iterative, sequences are one by one
     * aligned to cluster profile and added into the profile. If newId is
     * supported, result will be new cluster having this id. If nThreads is
     * supported, will run Clustal with --threads=nThreads option (not
     * recommended). Calls clusters setAsHasMSA() if successfull and returns
     * extended cluster. File with resulting MSA is saved in appropriate folder
     * (according to Settings) with name in format: [cluster ID].aln If a
     * cluster hasMSA() and a file with appropriate name already exists, new
     * alignment is constructed and old file is overwritten without warning.
     *
     *
     */
    @Override
    public ExtendClusterResult call() throws IOException, InterruptedException, DataException {
        List<UniqueSequence> addedSequences = new ArrayList<>();
        List<UniqueSequence> rejectedSequences = new ArrayList<>();
        ClustalRunner.multipleAlignment(extendedCluster);  //ensure cluster has MSA
        int maxAlnLength = Hammock.maxAlnLength;
        if (!Hammock.extensionIncreaseLength) {
            maxAlnLength = FileIOManager.getAlnLength(Settings.getInstance().getMsaDirectory() + extendedCluster.getId() + ".aln");
        }

        Cluster newCluster;
        newCluster = extendedCluster;
        for (int i = 0; i < insertedSequences.size(); i++) {

            if (!Statistics.checkCorrelation(newCluster, insertedSequences.get(i), Hammock.minCorrelation)) { //check correlation first
                rejectedSequences.add(insertedSequences.get(i));
            } else {

                String fastaPath = Settings.getInstance().getFastaDirectory() + newCluster.getId() + "_new.fa";
                FileIOManager.saveUniqueSequenceToFastaWithoutLabels(insertedSequences.get(i), fastaPath, newCluster.getId() + "_" + (extendedCluster.getUniqueSize() + i + 1));
                List<String> parameters = new ArrayList<>();
                parameters.add("--profile2");
                parameters.add(fastaPath);
                parameters.add("--profile1");
                parameters.add(Settings.getInstance().getMsaDirectory() + newCluster.getId() + ".aln");
                parameters.add("-o");
                parameters.add(Settings.getInstance().getMsaDirectory() + newCluster.getId() + "_testing.aln");
                parameters.add("--is-profile");
                parameters.add("--force");
                List<String> clustalParameters = Settings.getInstance().getClustalParameters();
                if (clustalParameters != null) {
                    parameters.addAll(clustalParameters);
                }
                ExternalProcessRunner.runProcess(Settings.getInstance().getClustalCommand(), parameters, System.out, System.err, "while extending cluster " + extendedCluster.getId() + " ");
                List<String> newClusterLines = FileIOManager.getAlignmentLines(Settings.getInstance().getMsaDirectory() + newCluster.getId() + "_testing.aln");
                if (FileIOManager.checkAlnLength(newClusterLines, maxAlnLength)
                        && FileIOManager.checkBothInnerGaps(newClusterLines, Hammock.maxInnerGaps)
                        && FileIOManager.checkConservedStates(newClusterLines, minMatchStates, minIc, Hammock.maxGapProportion)) {
                    FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + newCluster.getId() + "_testing.aln", Settings.getInstance().getMsaDirectory() + newCluster.getId() + ".aln");
                    addedSequences.add(insertedSequences.get(i));
                } else {
                    rejectedSequences.add(insertedSequences.get(i));
                }
            }
        }
        newCluster.insertAll(addedSequences);
        newCluster.setAsHasMSA();
        return new ExtendClusterResult(newCluster, rejectedSequences);
    }
}

/**
 * Class represents result from extendCluster routine
 *
 * @author Adam Krejci
 */
class ExtendClusterResult {

    private final Cluster cluster;
    private final List<UniqueSequence> rejectedSequences;

    public ExtendClusterResult(Cluster cluster, List<UniqueSequence> rejectedSequences) {
        this.cluster = cluster;
        this.rejectedSequences = rejectedSequences;
    }

    public Cluster getCluster() {
        return cluster;
    }

    public List<UniqueSequence> getRejectedSequences() {
        return rejectedSequences;
    }
}
