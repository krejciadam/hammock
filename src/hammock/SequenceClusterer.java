/*
 * Sequence clusterers are able to turn a collection of sequences into a a colleciton
   of clusters. 
 */

package hammock;

import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Adam Krejci
 */
public interface SequenceClusterer {
    /**
     * Clusters a collection of sequences using a sequence-sequence scoring scheme of choice.
     * @param sequences Sequences to be clustered.
     * @param scorer Scorer to be used to calculate sequence-sequence scores.
     * @return A collection of clusters.
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws DataException 
     */
    public List<Cluster> cluster(List<UniqueSequence> sequences, AligningSequenceScorer scorer) throws InterruptedException, ExecutionException, DataException;
    
}
