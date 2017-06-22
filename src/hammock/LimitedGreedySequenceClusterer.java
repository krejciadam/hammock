/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hammock;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author akrejci
 */
public class LimitedGreedySequenceClusterer implements SequenceClusterer {
    private final int threshold;
    private final int maxClusters;

    public LimitedGreedySequenceClusterer(int threshold, int maxClusters) {
        this.threshold = threshold;
        this.maxClusters = maxClusters;
    }
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences, AligningSequenceScorer sequenceScorer) throws InterruptedException, ExecutionException, DataException{
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
