/*
 * Performs cluster + database clustering using hmm-based clustering.
 */
package hammock;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 *
 * @author Adam Krejci
 */
public class IterativeHmmClusterer {

    /**
     * Performs iterative clustering using hmms. In each iteration, database
     * sequences are searched using cluster cores. Cores are then extended with
     * sequences found with score betted than actual value of
     * hmmsearchAssignThresholdSequence. Clusters having any overlap in theirs
     * sets of sequences found with score higher than hmmsearchOverlapThreshold
     * are then compared using hmm vs. hmm comparison and merged, if the score
     * of this comparison extends hhClusteringThreshold. Resulting clusters are
     * used as cluster cores for next iteration.
     *
     * @param clusterCores initial clusters to be extended
     * @param databaseSequences sequences to be iteratively searched and added
     * to clusters
     * @param hmmsearchAssignThresholdSequence array of minimal scores needed
     * for a sequence to be inserted into a cluster. Length of this array
     * determines the number of algorithm iterations.
     * @param hmmsearchOverlapThresholdSequence minimal score of a hmmsearch
     * sequence hit for appropriate clusters to be considered having overlapping
     * sets of sequences identified using hmmmsearch
     * @param hhClusteringThresholdSequence minimal hmm vs. hmm comparison score
     * needed for two clusters to be merged.
     * @param fullHHClustering sequence of booleans of length equal to the
     * number of rounds performed. When true, full HH clustering will be
     * performed in appropriate round
     * @param minMatchStates Minimal number of match states of a HMM after
     * merging
     * @param minMatchStatesSatisfyingIc Minimal number of match states having
     * information content higher than minIc
     * @param minIc minimal information content
     * @param nThreads Number of computational threads
     * @return List of resulting clusters
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws Exception
     */
    public static AssignmentResult iterativeHmmClustering(
            Collection<Cluster> clusterCores,
            Collection<UniqueSequence> databaseSequences,
            double[] hmmsearchAssignThresholdSequence,
            double[] hmmsearchOverlapThresholdSequence,
            double[] hhClusteringThresholdSequence,
            boolean[] fullHHClustering,
            Integer minMatchStates,
            Double minIc,
            int maxAlnLength,
            int nThreads) throws IOException, InterruptedException, ExecutionException, Exception {

        AssignmentResult currentState = new AssignmentResult(clusterCores, databaseSequences);
        for (int round = 0; round < hmmsearchAssignThresholdSequence.length; round++) {
            double hmmsearchAssignThreshold = hmmsearchAssignThresholdSequence[round];
            double hmmsearchOverlapThreshold = hmmsearchOverlapThresholdSequence[round];
            double hhClusteringThreshold = hhClusteringThresholdSequence[round];
            Hammock.logger.logAndStderr("");
            Hammock.logger.logAndStderr("Round " + (round + 1) + ":" + "\n");
            Hammock.logger.logAndStderr(currentState.getClusters().size() + " clusters remaining");
            Hammock.logger.logAndStderr("Building hmms and searching database...");
            List<HmmsearchSequenceHit> hits = HmmerRunner.searchWithHmms(currentState.getClusters(), currentState.getDatabaseSequences());
            Set<UnorderedPair<Cluster>> overlapingPairs = getOverlapingPairs(hits, hmmsearchOverlapThreshold); //must be before assignToClusters, because assignToClusters changes clusters
            Hammock.logger.logAndStderr("Extending clusters...");
            currentState = assignToClusters(currentState.getClusters(), hits, currentState.getDatabaseSequences(), hmmsearchAssignThreshold); //changes clusters
            Set<Set<Cluster>> mergeGroups;

            Set<Cluster> currentClusters = new HashSet<>();
            if (fullHHClustering[round] == false) {
                mergeGroups = IterativeHmmClusterer.getMergeGroups(overlapingPairs);
                currentClusters = new HashSet<>(currentState.getClusters());
                for (UnorderedPair<Cluster> pair : overlapingPairs) {
                    currentClusters.remove(pair.getBigger());
                    currentClusters.remove(pair.getSmaller());
                }
                Hammock.logger.logAndStderr(overlapingPairs.size() + " cluster pairs to check and merge.");
                Hammock.logger.logAndStderr("Merging clusters from " + mergeGroups.size() + " groups...");
            } else {
                currentClusters = new HashSet<>();
                mergeGroups = new HashSet<>();
                mergeGroups.add(new HashSet<>(currentState.getClusters()));
                Hammock.logger.logAndStderr("Overlap threshold is 0. Running full cluster merging.");
            }

            String path;
            File dir;
            if (!Hammock.inGalaxy) {
                path = Hammock.workingDirectory + Hammock.separatorChar + "alignments_other" + Hammock.separatorChar + "round_" + (round + 1) + "_after_assignment";
                dir = new File(path);
                dir.mkdir();
                for (Cluster cl : currentState.getClusters()) {
                    FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln", path + "/" + cl.getId() + ".aln");
                }
            }

            /*Build hhs in a well-paralelized fashion (before calling parallelHHClustering)*/
            Set<Cluster> clustersToMerge = new HashSet<>();
            for (Set<Cluster> mergeGroup : mergeGroups) {
                clustersToMerge.addAll(mergeGroup);
            }

            Hammock.logger.logAndStderr("Buiding hhs...");

            HHsuiteRunner.buildHHs(clustersToMerge, Hammock.threadPool);

            Hammock.logger.logAndStderr("HH clustering...");

            List<List<Cluster>> clusterLists = (parallelHHClustering(mergeGroups, hhClusteringThreshold, minMatchStates, minIc, maxAlnLength));
            for (List<Cluster> list : clusterLists) {
                currentClusters.addAll(list);
            }
            currentState = new AssignmentResult(new ArrayList<>(currentClusters), currentState.getDatabaseSequences());

            if (!Hammock.inGalaxy) {
                path = Hammock.workingDirectory + Hammock.separatorChar + "alignments_other" + Hammock.separatorChar + "round_" + (round + 1) + "_after_merging";
                dir = new File(path);
                dir.mkdir();
                for (Cluster cl : currentState.getClusters()) {
                    FileIOManager.copyFile(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln", path + "/" + cl.getId() + ".aln");
                }
            }

        }
        return currentState;
    }

    /**
     * Takes hmmsearch results as input. Inserts every sequence with any hits
     * having score >= scoreThreshold into a cluster with the highest score. If
     * two or more clusters have the same score, random cluster is selected out
     * of these. Returns resulting (extended) clusters and remaining sequences
     * (those not inserted into any cluster)
     *
     * @param clusters Clusters to be extended
     * @param hits Hmmsearch results
     * @param databaseSequences sequences to be inserted into clusters
     * @param scoreThreshold minimal hmmsearch score needed for a sequence to be
     * inserted into a cluster
     * @param nThreads Number of computational threads
     * @return AssignResult containing extended input clusters and sequences
     * that were not inserted into any cluster.
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     */
    private static AssignmentResult assignToClusters(Collection<Cluster> clusters, Collection<HmmsearchSequenceHit> hits, Collection<UniqueSequence> databaseSequences, Double scoreThreshold) throws IOException, InterruptedException, ExecutionException {
        Set<UniqueSequence> remainingSequences = new HashSet<>(databaseSequences);
        Set<Cluster> newClusters = new HashSet<>(clusters);

        /*Best hits finding:*/
        Map<UniqueSequence, HmmsearchSequenceHit> hitMap = new HashMap<>();
        for (HmmsearchSequenceHit hit : hits) {
            if (hitMap.containsKey(hit.getSequence())) {
                if (hitMap.get(hit.getSequence()).compareTo(hit) < 0) {
                    hitMap.put(hit.getSequence(), hit);
                }
            } else {
                if (hit.getScore() >= scoreThreshold) {
                    hitMap.put(hit.getSequence(), hit);
                }
            }
        }

        Hammock.logger.logAndStderr(hitMap.size() + " sequences to be inserted into clusters");

        /*Cluster extension: */
        Map<Cluster, List<HmmsearchSequenceHit>> extensionMap = new HashMap<>();
        for (Map.Entry<UniqueSequence, HmmsearchSequenceHit> entry : hitMap.entrySet()) {
            Cluster extendedCluster = entry.getValue().getCluster();
            List<HmmsearchSequenceHit> insertedSequences = extensionMap.get(extendedCluster);
            if ((insertedSequences) == null) {
                insertedSequences = new ArrayList<>();
            }
            insertedSequences.add(entry.getValue());
            extensionMap.put(extendedCluster, insertedSequences);
            remainingSequences.remove(entry.getKey());
        }

        Hammock.logger.logAndStderr(extensionMap.size() + " clusters to be extended");

        ExtendedClusters extendedClusters = ClustalRunner.extendClusters(extensionMap, Hammock.minMatchStates, Hammock.minIc);
        newClusters.addAll(extendedClusters.getClusters());
        newClusters.addAll(clusters);
        Hammock.logger.logAndStderr(extendedClusters.getRejectedSequences().size() + " sequences rejected");
        remainingSequences.addAll(extendedClusters.getRejectedSequences());
        return new AssignmentResult(newClusters, remainingSequences);
    }

    /**
     * Returns all pairs of clusters that have any overlap in their sets of
     * sequences found by hmmsearch. i.e. for any two clusters A and B, if there
     * is at least one sequence S such that score(A, S) >= scoreThreshold and
     * score(B, S) >= threshold, UnorderedPair containing A and B is reported.
     *
     * @param hits results of hmmsearch to be processed
     * @param scoreThreshold minimal score for sequences to be considered as
     * pair of the overlap.
     * @return All pairs of Clusters having any overlap in their sets of
     * sequences found by hmmsearch.
     */
    private static Set<UnorderedPair<Cluster>> getOverlapingPairs(Collection<HmmsearchSequenceHit> hits, Double scoreThreshold) {
        Map<UniqueSequence, List<Cluster>> hitMap = new HashMap<>();
        for (HmmsearchSequenceHit hit : hits) {
            if (hit.getScore() >= scoreThreshold) {
                List<Cluster> newList;
                if (hitMap.containsKey(hit.getSequence())) {
                    newList = hitMap.get(hit.getSequence());
                } else {
                    newList = new ArrayList<>();
                }
                newList.add(hit.getCluster());
                hitMap.put(hit.getSequence(), newList);
            }
        }
        Set<UnorderedPair<Cluster>> result = new HashSet<>();
        for (List<Cluster> clusterList : hitMap.values()) {
            if (clusterList.size() >= 2) {
                for (int i = 0; i < clusterList.size() - 1; i++) {
                    for (int j = i + 1; j < clusterList.size(); j++) {
                        result.add(new UnorderedPair<>(clusterList.get(i), clusterList.get(j)));
                    }
                }
            }
        }
        return result;
    }

    /**
     * Parses graph of cluster pairs into set of (weakly) connected components
     * using DFS
     *
     * @param pairs Pairs of clusters (== list of graph edges, in the graph
     * meaning).
     * @return Set of groups of clusters. Each group represents one (weakly)
     * connected component
     */
    private static Set<Set<Cluster>> getMergeGroups(Collection<UnorderedPair<Cluster>> pairs) {
        Map<Cluster, Set<Cluster>> nodeMap = new HashMap<>();
        for (UnorderedPair<Cluster> pair : pairs) {
            Set<Cluster> smallerSet = nodeMap.get(pair.getSmaller());
            Set<Cluster> biggerSet = nodeMap.get(pair.getBigger());
            if (smallerSet == null) {
                smallerSet = new HashSet<>();
            }
            if (biggerSet == null) {
                biggerSet = new HashSet<>();
            }
            smallerSet.add(pair.getBigger());
            biggerSet.add(pair.getSmaller());
            nodeMap.put(pair.getSmaller(), smallerSet);
            nodeMap.put(pair.getBigger(), biggerSet);
        }
        Set<Set<Cluster>> result = new HashSet<>();
        while (nodeMap.size() > 0) {
            Cluster lineCluster = nodeMap.keySet().iterator().next();
            result.add(IterativeHmmClusterer.processMapLine(nodeMap, lineCluster));
        }
        return result;
    }

    /**
     * Recursive method for DFS
     *
     * @param nodeMap
     * @param lineCluster
     * @return
     */
    private static Set<Cluster> processMapLine(Map<Cluster, Set<Cluster>> nodeMap, Cluster lineCluster) {
        Set<Cluster> result = new HashSet<>();
        Set<Cluster> toProcess = nodeMap.get(lineCluster);
        if (toProcess == null) {
            return result;
        }
        nodeMap.remove(lineCluster);
        for (Cluster cl : toProcess) {
            result.add(cl);
            result.addAll(IterativeHmmClusterer.processMapLine(nodeMap, cl));
        }
        return result;
    }

    /**
     * Runs hhClustering in multiple instances (one for each clusterGroup), i.e.
     * performs hh clustering within each clusterGroup. Number of
     * computation-intensive threads will never exceed nThreads.
     *
     * @param clusterGroups groups of clusters to be clusterer.
     * @param scoreThreshold Minimal hh alignment score for 2 clusters to be
     * merged
     * @param nThreads Number of computational threads
     * @return
     * @throws InterruptedException
     * @throws ExecutionException
     */
    private static List<List<Cluster>> parallelHHClustering(Set<Set<Cluster>> clusterGroups, Double scoreThreshold, Integer minMatchStates, Double minIc, int maxAlnLength) throws InterruptedException, ExecutionException {
        ExecutorService lowThreadPool = Executors.newFixedThreadPool(Hammock.nThreads);
        CompletionService<List<Cluster>> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);

        for (Collection<Cluster> clusterGroup : clusterGroups) {
            resultPool.submit(new hhClusteringRunner(clusterGroup, scoreThreshold, minMatchStates, minIc, maxAlnLength, lowThreadPool));
        }
        List<List<Cluster>> result = new ArrayList<>();
        for (int i = 0; i < clusterGroups.size(); i++) {
            result.add(resultPool.take().get());
        }
        lowThreadPool.shutdown();
        return result;
    }

    /**
     * Perform clustering based on hmm-hmm comparison. At each round, 2 most
     * similar clusters are merged. Merging continues until there are no
     * clusters having alignment score >= scoreThreshold.
     *
     * @param initialClusters Clusters to be clustered
     * @param scoreThreshold Minimal hh alignment score for 2 clusters to be
     * merged
     * @param minMatchStates Minimal number of HMM match states after hmm
     * merging
     * @param lowThreadPool an ExecutorService to be used for whole process of
     * clustering. May be shared between several instances of hhClustering
     * @param mergeThreads
     * @return Resulting List of (merged) clusters
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws IOException
     * @throws DataException
     */
    public static List<Cluster> hhClustering(Collection<Cluster> initialClusters, Double scoreThreshold, Integer minMatchStates, Double minIc, int maxAlnLength, ExecutorService lowThreadPool) throws InterruptedException, ExecutionException, IOException, DataException, Exception {
        List<Cluster> clusterList = new ArrayList<>();
        clusterList.addAll(initialClusters);

        boolean build = false;  //chceck if all clusters have their HHs constructed (they  shoud...)
        for (Cluster cl : initialClusters) {
            if (cl.hasHH()) {
                build = true;
                break;
            }
        }
        if (build) {
            HHsuiteRunner.buildHHs(clusterList, lowThreadPool); //build all hhs here (building from multiple threads in next step is not desired)
        }
        /*initial all-vs-all comparison*/
        List<HHalignHit> toAdd = new ArrayList<>();
        for (HHalignHit hit : HHsuiteRunner.alignAllVsAll(clusterList, lowThreadPool)) {
            if (hit.getScore() >= scoreThreshold) {
                toAdd.add(hit);
            }
        }

        SortedSet<HHalignHit> hitSet;
        if (toAdd.size() > 0) { //we continue only if there is at least one good hit
            hitSet = new TreeSet<>(toAdd);
            HHalignHit firstPair = hitSet.last();
            while ((hitSet.size() >= 1) && (firstPair.getScore() >= scoreThreshold)) {
                hitSet.remove(firstPair); //remove the best scoring pair
                Cluster tempCluster = HHsuiteRunner.mergeClusters(firstPair, initialClusters.hashCode()); //temp cluster id = hash code. Needed for this thread to have unique temp cluster id
                HHsuiteRunner.builHH(tempCluster);
                if (Statistics.checkCorrelation(firstPair.getSearchedCluster(), firstPair.getFoundCluster(), Hammock.minCorrelation)
                        && (FileIOManager.checkMatchStatesAndIc(Settings.getInstance().getMsaDirectory() + tempCluster.getId() + ".aln", minMatchStates, minIc, Hammock.maxGapProportion))
                        && (FileIOManager.checkBothInnerGaps(Settings.getInstance().getMsaDirectory() + tempCluster.getId() + ".aln", Hammock.maxInnerGaps))
                        && (FileIOManager.checkAlnLength(Settings.getInstance().getMsaDirectory() + tempCluster.getId() + ".aln", maxAlnLength))) { //satisfies conditions
                    clusterList.remove(firstPair.getBiggerCluster());
                    clusterList.remove(firstPair.getSmallerCluster());
                    Set<HHalignHit> newHitSet = new HashSet<>();
                    for (HHalignHit pair : hitSet) { //remove all pairs containing any cluster from best scoring pair
                        if ((pair.getBiggerCluster().equals(firstPair.getBiggerCluster()))
                                || (pair.getBiggerCluster().equals(firstPair.getSmallerCluster()))
                                || (pair.getSmallerCluster().equals(firstPair.getBiggerCluster()))
                                || (pair.getSmallerCluster().equals(firstPair.getSmallerCluster()))) {
                        } else {
                            newHitSet.add(pair);
                        }
                    }
                    hitSet = new TreeSet<>(newHitSet);
                    Cluster newCluster = HHsuiteRunner.mergeClusters(firstPair, firstPair.getBiggerCluster().getId());
                    HHsuiteRunner.builHH(newCluster);
                    for (HHalignHit hit : HHsuiteRunner.alignHmmList(newCluster, clusterList, lowThreadPool)) {
                        if (hit.getScore() >= scoreThreshold) {
                            hitSet.add(hit);
                        }
                    }
                    clusterList.add(newCluster);
                }
                if (hitSet.size() > 0) {
                    firstPair = hitSet.last();
                }
            }
        }
        return clusterList;
    }

}

class hhClusteringRunner implements Callable<List<Cluster>> {

    private final Collection<Cluster> clusters;
    private final Double scoreThreshold;
    private final ExecutorService lowThreadPool;
    private final Integer minMatchStates;
    private final Double minIc;
    private final int maxAlnLength;

    public hhClusteringRunner(Collection<Cluster> clusters, Double scoreThreshold, Integer minMatchStates, Double minIc, int maxAlnLength, ExecutorService lowThreadPool) {
        this.clusters = clusters;
        this.scoreThreshold = scoreThreshold;
        this.lowThreadPool = lowThreadPool;
        this.minMatchStates = minMatchStates;
        this.minIc = minIc;
        this.maxAlnLength = maxAlnLength;

    }

    @Override
    public List<Cluster> call() throws Exception {
        return IterativeHmmClusterer.hhClustering(clusters, scoreThreshold, minMatchStates, minIc, maxAlnLength, lowThreadPool);
    }

}
