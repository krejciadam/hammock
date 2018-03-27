/*
 * Runs clinkage/greedy clustering within pre-defined groups
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
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
public class inGroupClusterer {
    
    private final int threshold;
    private final SequenceScorer scorer;
    private final CompletionService<List<Cluster>> resultPool;

    public inGroupClusterer(int threshold, SequenceScorer scorer,  ExecutorService threadPool) {
        this.threshold = threshold;
        this.scorer = scorer;
        resultPool = new ExecutorCompletionService<>(threadPool);
    }
    
    
    
    
    public List<Cluster> clusterWithinGroups(Collection<Cluster> groups) throws InterruptedException, ExecutionException{
        
        
        for (Cluster group : groups){
            SequenceClusterer clinkageClusterer = new ClinkageSequenceClusterer(
                scorer, threshold, Integer.MAX_VALUE, Executors.newSingleThreadExecutor()); //the max value means no caching will be performed
            resultPool.submit(new ClusteringRunner(clinkageClusterer, group.getSequences()));
        }
        List<Cluster> result = new ArrayList<>();
        for (Cluster group : groups){
            result.addAll(resultPool.take().get());
        }
        return(result);
        
    }
    
}

class ClusteringRunner implements Callable<List<Cluster>>{

    private final SequenceClusterer clusterer;
    private final List<UniqueSequence> sequences;

    public ClusteringRunner(SequenceClusterer clusterer, List<UniqueSequence> sequences) {
        this.clusterer = clusterer;
        this.sequences = sequences;
    }
    
    @Override
    public List<Cluster> call() throws Exception {
        List<Cluster> clusters = clusterer.cluster(sequences);
        return(clusters);
    }
    
}
