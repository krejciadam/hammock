/*
 * Clusters sequences using a greedy algorithm. The results satisfy the 
complete-likage criterion.
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class LimitedGreedySequenceClusterer implements SequenceClusterer {
    private final int threshold;
    private final int maxClusters;
    private final AligningSequenceScorer sequenceScorer;

    public LimitedGreedySequenceClusterer(AligningSequenceScorer sequenceScorer, int threshold, int maxClusters) {
        this.threshold = threshold;
        this.maxClusters = maxClusters;
        this.sequenceScorer = sequenceScorer;
    }
    
    /**
     * Clusters sequences using a greedy scheme. The results satisfy the complete-linkage
     * criterion.
     * @param sequences Sequences to cluster.
     * @return Resulting clustering. Orphan sequences are returned as single-member
     * clusters.
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws DataException 
     */
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences) throws InterruptedException, ExecutionException, DataException{
        ClusterScorer clusterScorer = new ClinkageClusterScorer(sequenceScorer, threshold);
        List<Cluster> clusters = firstPhase(sequences, maxClusters, clusterScorer);
        //now clusters (size > 1) are at the beginning of the list
        int index = clusters.size();
        for(int i = 0; i < clusters.size(); i++){
            if (clusters.get(i).getUniqueSize() == 1){
                index = i;
                break;
            }
        }
        List<Cluster> actualClusters = new ArrayList<>(clusters.subList(0, index));
        List<Cluster> actualSequences = new ArrayList<>(clusters.subList(index, clusters.size()));
        List<Cluster> remainingSequences = new ArrayList<>();
        CompletionService<NearestCluster> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        
        int sumCommodity = 0;
        for (Cluster cl : actualClusters){
            sumCommodity += cl.getUniqueSize();
        }
        for (Cluster cl : actualSequences){
            NearestCluster foundCluster = ClinkageSequenceClusterer.findNearestClusterParallel(actualClusters, cl, clusterScorer, sumCommodity, resultPool);
            if (foundCluster != null && foundCluster.getScore() >= threshold){
                foundCluster.getCluster().insertAll(cl.getSequences());
            }else{
                remainingSequences.add(cl);
            }
        }
        actualClusters.addAll(remainingSequences);
        return(actualClusters);
    }
    
    /**
     * First, sequences are merged with their nearest neighbors. Clusters of 
     * size more than 2 can arise, in case an existing cluster is a sequence's 
     * nearest neighbor. Stops when there are maxClusters clusters of size at least
     * 2 or all the sequences were processed.
     */
    private List<Cluster> firstPhase(List<UniqueSequence> sequences, int maxClusters, ClusterScorer clusterScorer) throws InterruptedException, ExecutionException, DataException{
        List<Cluster> initialList = new ArrayList<>();
        for (int i = 0; i < sequences.size(); i++){
            List<UniqueSequence> l = new ArrayList<>();
            l.add(sequences.get(i));
            initialList.add(new Cluster(l, i));
        }
        CompletionService<NearestCluster> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        
        List<Cluster> actualClusters = new ArrayList<>();
        List<Cluster> actualSequences = new ArrayList<>();
        int sumCommodityClusters = 0;
        int index = 0;
        while(index < initialList.size() && actualClusters.size() < maxClusters){
            Cluster comparedCluster = initialList.get(index);
            NearestCluster foundClusterClusters = ClinkageSequenceClusterer.findNearestClusterParallel(actualClusters, comparedCluster, clusterScorer, sumCommodityClusters, resultPool);
            NearestCluster foundClusterSequences = ClinkageSequenceClusterer.findNearestClusterParallel(initialList.subList(index + 1, initialList.size()), comparedCluster, clusterScorer, initialList.size() - index - 1, resultPool);
            if (foundClusterClusters != null){
                if(foundClusterSequences != null) {
                    if (foundClusterClusters.getScore() >= foundClusterSequences.getScore()){
                        foundClusterClusters.getCluster().insertAll(comparedCluster.getSequences());
                    } else{
                        comparedCluster.insertAll(foundClusterSequences.getCluster().getSequences());
                        actualClusters.add(comparedCluster);
                        initialList.remove(foundClusterSequences.getCluster());
                    }
                } else{
                    foundClusterClusters.getCluster().insertAll(comparedCluster.getSequences());
                }
            } else{
                if (foundClusterSequences != null){
                    comparedCluster.insertAll(foundClusterSequences.getCluster().getSequences());
                    actualClusters.add(comparedCluster);
                    initialList.remove(foundClusterSequences.getCluster());
                } else{
                    actualSequences.add(comparedCluster);
                }
            }
            index++;
        }
        actualClusters.addAll(actualSequences);
        actualClusters.addAll(initialList.subList(index, initialList.size()));
        return(actualClusters);
    }
}
