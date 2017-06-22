/*   Scorer returning seuquence score based on a scoring matrix. Sequences
 are aligned without any gaps, but shift relative to each other's start
 is allowed. A shift penalty is given for such shift. Score is the maximmum
 of all possible shifts. 
 */
package hammock;

/**
 *
 * @author Adam Krejci
 */
public class ShiftedScorer implements AligningSequenceScorer {

    private final int[][] scoringMatrix;
    private final int shiftPenalty;
    private final int maxShift;

    /**
     *
     * @param scoringMatrix a scoring matrix in correct format (see
     * FileIOManeger)
     * @param maxShift maximum number of positions the sequences will be shifted
     * by. shifts to both sides are performed.
     * @param shiftPenalty a (most likely negative) number to be added to score.
     * (shiftPenalty)*(actual shift) will be added to score examined for each
     * possible shift.
     */
    public ShiftedScorer(int[][] scoringMatrix, int shiftPenalty, int maxShift) {
        this.scoringMatrix = scoringMatrix;
        this.maxShift = maxShift;
        this.shiftPenalty = shiftPenalty;
    }

    /**
     * Returns alignment score based on a scoring matrix. Sequences are aligned
     * without any gaps, but may be shifted relative to each other's start. A
     * penalty is added for each AA aligned to a gap.
     *
     *
     * Returns the best hit along with second compared sequence.
     *
     * @param seq1 First sequence to be aligned
     * @param seq2 Second sequence to be aligned
     * @return
     * @throws hammock.DataException
     */
    @Override
    public AligningScorerResult scoreWithShift(UniqueSequence seq1, UniqueSequence seq2) throws DataException {
        UniqueSequence shorterSequence;
        UniqueSequence longerSequence;
        if (seq1.getSequence().length >= seq2.getSequence().length) {
            shorterSequence = seq2;
            longerSequence = seq1;
        } else {
            shorterSequence = seq1;
            longerSequence = seq2;
        }

        if (maxShift >= shorterSequence.getSequence().length) {
            throw new DataException("Shift too big: " + (shorterSequence.getSequence().length - 1)
                    + " is maximum, but " + maxShift + " found");
        }

        int bestScore = Integer.MIN_VALUE;
        int bestScoreShift = 0;
        int lengthDifference = (longerSequence.getSequence().length - shorterSequence.getSequence().length);
        for (int actualShift = -maxShift; actualShift <= maxShift + (lengthDifference); actualShift++) {
            int actualScore = 0;
            if (actualShift <= 0) { //first(shorter) sequence shifted to the left or not shifted
                for (int i = 0; i < shorterSequence.getSequence().length + actualShift; i++) {
                    actualScore += letterScore(shorterSequence.getSequence()[i - actualShift], longerSequence.getSequence()[i]);
                }
            } else { //first(shorter) sequence shifted to the right
                for (int i = 0; i < Math.min(shorterSequence.getSequence().length, longerSequence.getSequence().length - actualShift); i++) {
                    actualScore += letterScore(shorterSequence.getSequence()[i], longerSequence.getSequence()[i + actualShift]);
                }
            }
            //actualScore += Math.abs(actualShift) * shiftPenalty;
            actualScore += (lengthDifference) * shiftPenalty; //these gaps are always present
            if (actualShift < 0){
                actualScore += -actualShift * 2 * shiftPenalty; //shorter shifted to the left
            }
            if (actualShift > lengthDifference){
                actualScore += (actualShift - lengthDifference) * 2 * shiftPenalty; //shorter shifted to the right
            }
            if (actualScore > bestScore) {
                bestScore = actualScore;
                bestScoreShift = actualShift;
            }
        }
        if (shorterSequence != seq2){
            bestScoreShift = -bestScoreShift; //best shift was calculated as shift of shorter sequence
        }
        return new AligningScorerResult(bestScore, bestScoreShift, seq2);
    }

    @Override
    public int sequenceScore(UniqueSequence seq1, UniqueSequence seq2) throws DataException {
        return scoreWithShift(seq1, seq2).getScore();
    }

    /**
     * Returns letter to letter score based on scoring matrix provided
     *
     * @param firstLetter one (any) of letters to be compared
     * @param secondLetter one (any) of letter to be compared
     * @param scoringMatrix scoring matrix in correct format (see FileIOManager)
     * @return
     */
    private int letterScore(int firstLetter, int secondLetter) {
        return scoringMatrix[firstLetter][secondLetter];
    }

}
