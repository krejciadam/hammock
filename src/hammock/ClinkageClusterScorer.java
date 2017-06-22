/*
 * Performs cluster-cluster comparison
 */
package hammock;

import java.util.Collection;

/**
 *
 * @author Adam Krejci
 */
public class ClinkageClusterScorer implements ClusterScorer {
    private final SequenceScorer scorer;
    private final int threshold;

    public ClinkageClusterScorer(SequenceScorer scorer, int threshold) {
        this.scorer = scorer;
        this.threshold = threshold;
    }
    
    /**
     * Returns complete-linkage score between 2 clusters or Integer.min_value 
     * when the score is less than threshold
     * @param cl1 First (any) cluster to be compared
     * @param cl2 Second (any) cluster to be compared
     * @param scorer Scorer to be used for sequence-sequence comparison
     * @param threshold If cluster-cluster score is less than this threshold, 
     * it will not be calculated and reported. Integer.min_value will be reported instead
     * @return Complete-linkage score of the clusters
     * @throws hammock.DataException
     */
    @Override
    public int clusterScore(Cluster cl1, Cluster cl2) throws DataException{
        Collection<UniqueSequence> sequences1 = cl1.getSequences();
        Collection<UniqueSequence> sequences2 = cl2.getSequences();
        
        int result = Integer.MAX_VALUE;
        
        for (UniqueSequence seq1 : sequences1){
            for (UniqueSequence seq2 : sequences2){
                int roundRes = scorer.sequenceScore(seq1, seq2);
                if (roundRes < result){
                    result = roundRes;
                    if (result < threshold){
                        return(Integer.MIN_VALUE + 1);
                    }
                }
                
            }
        }
        return(result);
    }
    
}