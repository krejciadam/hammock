/*
 * Class carries the result of shifted score - actual score and shift used to achieve it
 */

package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class AligningScorerResult implements Comparable<AligningScorerResult> {
    
    private final int score;
    private final int shift;
    private final UniqueSequence sequence;

    public AligningScorerResult(int score, int shift, UniqueSequence sequence) {
        this.score = score;
        this.shift = shift;
        this.sequence = sequence;
    }

    public int getScore() {
        return score;
    }

    public int getShift() {
        return shift;
    }

    public UniqueSequence getSequence() {
        return sequence;
    }

    @Override
    public int compareTo(AligningScorerResult o) {
        if (this.getScore() == o.getScore()){
            return this.getSequence().compareTo(o.getSequence());
        } else{
            return this.getScore() - o.getScore();
        }
    }
    
}
    
