/*
 * Sequence clusterers are able to turn a collection of sequences into a a colleciton
   of clusters. 
 */

package hammock;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Adam Krejci
 */
public interface SequenceClusterer {
    
    public List<Cluster> cluster(Collection<UniqueSequence> sequences, ShiftedScorer scorer) throws InterruptedException, ExecutionException;
    
}
