/*
 * Interface for sequence scorers perforimg sequence alignment
 */
package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public interface AligningSequenceScorer extends SequenceScorer {
     public AligningScorerResult scoreWithShift(UniqueSequence seq1, UniqueSequence seq2)  throws DataException;
}
