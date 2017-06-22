/*
 * Two sequences may be compared in several ways. Scorer objects take two sequences
   and return their similarity as an integer score.
 */

package hammock;

/**
 *
 * @author Adam Krejci
 */
public interface SequenceScorer {
    
    public int sequenceScore(UniqueSequence seq1, UniqueSequence seq2) throws DataException;
}
