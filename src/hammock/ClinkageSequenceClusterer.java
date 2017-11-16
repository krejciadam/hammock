/*
 * Clusters sequences using the complete-linkage scheme.
 */
package hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class ClinkageSequenceClusterer implements SequenceClusterer {


    private final int sizeLimit; //only put clusters of this size and larger into cache.
                                //(only applies for CachedCLScorer)
    private final int threshold;

    public ClinkageSequenceClusterer(int threshold) {
        this.threshold = threshold;
        this.sizeLimit = 1;
    }
    
    public ClinkageSequenceClusterer(int threshold, int sizeLimit) {
        this.threshold = threshold;
        this.sizeLimit = sizeLimit;
    }
     
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences, AligningSequenceScorer sequenceScorer) throws InterruptedException, ExecutionException {
        CachedClusterScorer clusterScorer = new CachedClusterScorer(new ClinkageClusterScorer(sequenceScorer, threshold), sizeLimit);
        ClusterStack<Cluster> stack = new ClusterStack<>();
        Set<Cluster> readyClusters = new HashSet<>();  //a list for clusters whose nearest neighbor is past the threshold
        int sumCommodity = 0; //Sum of price of all the remaining clusters    
        Set<Cluster> activeClusters = new HashSet<>();
        int currentId = 1;
        for (UniqueSequence seq : sequences){
            Set<UniqueSequence> seqs = new HashSet<>();
            seqs.add(seq);
            activeClusters.add(new Cluster(seqs, currentId));
            currentId ++;
        }
        
        for (Cluster cluster : activeClusters) {
            sumCommodity += cluster.getUniqueSize();
        }
        
        CompletionService<NearestCluster> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        
        while (activeClusters.size() > 1) {
            if (!Hammock.inGalaxy){
                System.err.print("\rSequences remaining to process: "
                        + activeClusters.size() + ", clusters ready: "
                        + readyClusters.size() + " ");
            }

            Cluster randomCluster = activeClusters.iterator().next();
            stack.push(randomCluster); // choose arbitrary cluster
            while (!stack.isEmpty()) {
                Cluster top = stack.peek();
                
                /*Nearest neighbor finding:*/
                
                NearestCluster foundNearestCluster = findNearestClusterParallel(activeClusters, top, clusterScorer, sumCommodity, resultPool);
                int maxScore = Integer.MIN_VALUE;
                Cluster nearestCluster = null;
                if (foundNearestCluster != null){
                    nearestCluster = foundNearestCluster.getCluster();
                    maxScore = foundNearestCluster.getScore();
                }
                
                 /*If the score is too low, remove top and start over (this can only happen when there is only one cluster in the stack)*/
                if (maxScore < threshold) {
                    stack.pop();
                    readyClusters.add(top);
                    activeClusters.remove(top);
                    sumCommodity -= top.getUniqueSize();
                    continue;
                }
                
                /*The nearest cluster to "top" is the 2nd in the stack
                 * Remove both from the stack, merge them and put the new into active clusters*/
                if ((stack.size() > 1) && stack.peekSecond().equals(nearestCluster)) { 
                    currentId++;
                    stack.pop();
                    stack.pop();
                    activeClusters.remove(top);
                    activeClusters.remove(nearestCluster);
                    clusterScorer.join(top, nearestCluster, currentId);
                    sumCommodity -= (top.getUniqueSize() + nearestCluster.getUniqueSize()); 
                    
                    List<UniqueSequence> mergedSeqs = top.getSequences();
                    mergedSeqs.addAll(nearestCluster.getSequences());
                    
                    Cluster newTop = new Cluster(mergedSeqs, currentId);

                    activeClusters.add(newTop);
                    sumCommodity += newTop.getUniqueSize();
                } else { //only 1 clyuster in the stack or 2nd is not equal to the closest
                    stack.push(nearestCluster);
                }
            }
        }
        /*The last cluster is also ready*/
        readyClusters.add(activeClusters.iterator().next());

        /*return list*/
        List<Cluster> toReturn = new ArrayList<>();
        toReturn.addAll(readyClusters);
        return new ArrayList<>(toReturn);
    }
    
    /**
     * Finds the most similar cluster in a collection. 
     * @param inputClusters the collection to be searched
     * @param comparedCluster the cluster to find the nearest neighbor of
     * @param scorer Scoring scheme to use for cluster-cluster comparison
     * @param sumCommodity Actual sum of cluster calculation price
     * @param resultPool To be used for parallelization
     * @return The best scoring cluster and comparison information.
     * @throws InterruptedException
     * @throws ExecutionException 
     */
    public static NearestCluster findNearestClusterParallel(Collection<Cluster> inputClusters, Cluster comparedCluster, ClusterScorer scorer, int sumCommodity, CompletionService<NearestCluster> resultPool) throws InterruptedException, ExecutionException{
        if(inputClusters.isEmpty()){
            return(new NearestCluster(null, Integer.MIN_VALUE));
        }
        int currentMaxNumberOfThreads = ClinkageSequenceClusterer.setNumberOfParts(inputClusters);
        List<Set<Cluster>> activeClustersParts = ClinkageSequenceClusterer.activeClustersParts(inputClusters, sumCommodity, currentMaxNumberOfThreads);
        /*Vlastni hledani nejblizsiho clusteru*/
        for (int i = 0; i < activeClustersParts.size(); i++) {
            resultPool.submit(
                new NearestClusterRunner(
                activeClustersParts.get(i),
                comparedCluster, scorer));
            }
                
        int maxScore = (Integer.MIN_VALUE + 42); //If NearestCluster from the next part is NULL, it has MIN_VALUE. This way, we find out that these cases are skipped right away.
        NearestCluster nearestCluster = null;
                
            /*Results from the threads that already finished*/
        for (Set<Cluster> activeClustersPart : activeClustersParts) {
            NearestCluster currentCluster = resultPool.take().get();     
            
            
            if (currentCluster.getScore() < maxScore) {
                continue;
            }

            if (currentCluster.getScore() > maxScore) {
                nearestCluster = currentCluster;
                maxScore = nearestCluster.getScore();
            } else { //currentScore == maxScore   
                if ((currentCluster.getCluster().size() > nearestCluster.getCluster().size())) { //when currentCluster is null (exceptional situation. Onlu the top was in the thread), we can skip
                    nearestCluster = currentCluster;
                } else{
                    if ((currentCluster.getCluster().size() == nearestCluster.getCluster().size()) && (currentCluster.getCluster().getId() < nearestCluster.getCluster().getId())){ //equal size
                        nearestCluster = currentCluster;
                    }
                }
            }
        }
        return(nearestCluster);
    }
    
    
     /**
     * Finds how many parts to divide the work into
     *
     * @param activeClusters 
     * @return
     */
    private static int setNumberOfParts(Collection<Cluster> activeClusters) {
        int currentMaxNumberOfParts = Hammock.nThreads * 4;
        if (activeClusters.size() < Hammock.nThreads * 4 + 1) {      
            currentMaxNumberOfParts = Math.max(activeClusters.size() - 1, 1);
        }
        return currentMaxNumberOfParts;
    }
    
     /**
     * Divide the work between threads. Ideally evenly.
     *
     * @param activeClusters A collection to divide in the parts
     * @param sumCommodity Sum of "price" (i.e. expected calculation time)
     * @param numberOfParts Divide into this many parts
     * @return Each set in the list is a portion for 1 thread.
     */
    private static List<Set<Cluster>> activeClustersParts(Collection<Cluster> activeClusters,
            long sumCommodity,
            int numberOfParts) {

        List<Set<Cluster>> result = new ArrayList<>();
        Set<Cluster> currentSet = new HashSet<>();
        long forOneThread = (sumCommodity / numberOfParts) + 1; //DIV
        long portionCounter = forOneThread;
        for (Cluster clust : activeClusters) {
            currentSet.add(clust);
            portionCounter -= clust.getUniqueSize();
            if (portionCounter <= 0) {
                result.add(currentSet);
                currentSet = new HashSet<>();
                portionCounter = forOneThread;
            }
        }
        if (currentSet.size() > 0) {
            result.add(currentSet); 
        }
        return result;
    }
    
    
}

class ClusterStack<E> extends Stack<E> {
    /**
     * Returns the second element in the stack
     * @return The second stack element
     */
    public E peekSecond() {
        E element = this.get(this.size() - 2);
        return element;
    }
}


/**
 * For parallel search for nearest clusters
 */
class NearestClusterRunner implements Callable<NearestCluster> {

    private final Collection<Cluster> databaseClusters;
    private final Cluster comparedCluster;
    private final ClusterScorer scorer;

    public NearestClusterRunner(Collection<Cluster> activeClustersPart,
            Cluster top,
            ClusterScorer scorer) {
        this.databaseClusters = activeClustersPart;
        this.comparedCluster = top;
        this.scorer = scorer;
    }

    @Override
    public NearestCluster call() throws DataException{
        int maxScore = Integer.MIN_VALUE;
        int score;
        Cluster nearestCluster = null;
        for (Cluster i : databaseClusters) {
            score = scorer.clusterScore(i, comparedCluster);
            if (score < maxScore) {
                continue;
            }
            if (score > maxScore) {
                if ((i != comparedCluster)) {
                    maxScore = score;
                    nearestCluster = i;
                }
            } /*else: (score == maxScore) If there are more most similar clusters,      SHOULD IMPEMENT A COMPARATOR HERE
            we need a DETERMINISTIC decision which is the most similar. We use size
            and id. 
             */ else {
                if (i != comparedCluster){
                    if (i.size() > nearestCluster.size()){
                        nearestCluster = i;                        
                    } else{
                        if (i.size() < nearestCluster.size()){
                           //nothing, no change 
                        }
                        else { //size is the same as well. Decide using id.
                            if ((i.getId() < nearestCluster.getId())) {
                            nearestCluster = i;
                        }    
                        }
                    }
                }
            }
        }
        return new NearestCluster(nearestCluster, maxScore);
    }
}