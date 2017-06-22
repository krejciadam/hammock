/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hammock;

/**
 *
 * @author akrejci
 */
public interface AligningSequenceScorer extends SequenceScorer {
     public AligningScorerResult scoreWithShift(UniqueSequence seq1, UniqueSequence seq2)  throws DataException;
}
