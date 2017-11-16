/*
 * Interface for all clases calculating cluster-cluster score
 */
package hammock;

/**
 *
 * @author Adam Krejci
 */
public interface ClusterScorer {
    /**
     * Returns the score of two clusters
     * @param cl1 one (any) cluster to be scored
     * @param cl2 one (any) cluster to be scored
     * @return cluster-cluster score.
     * @throws DataException 
     */
    public int clusterScore(Cluster cl1, Cluster cl2) throws DataException;
    
}
