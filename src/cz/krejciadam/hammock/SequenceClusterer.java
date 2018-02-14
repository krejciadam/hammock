/*
 * Sequence clusterers are able to turn a collection of sequences into a a colleciton
   of clusters. 
 */

package cz.krejciadam.hammock;

import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Adam Krejci
 */
public interface SequenceClusterer {
    /**
     * Clusters a collection of sequences.
     * @param sequences Sequences to be clustered.
     * @return A collection of clusters.
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws DataException 
     */
    public List<Cluster> cluster(List<UniqueSequence> sequences) throws InterruptedException, ExecutionException, DataException;
    
}
