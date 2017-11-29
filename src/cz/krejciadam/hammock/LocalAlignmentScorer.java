/*
 * This class performs local alignment using the Smith-Waterman algorithm.
 */
package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class LocalAlignmentScorer implements SequenceScorer {

    private final int[][] scoringMatrix;
    private final int gapOpenPenalty;
    private final int gapExtendPenalty;

    private enum Direction {
        LEFT, UP, DIAGONAL, NOWHERE;
    }

    public LocalAlignmentScorer(int[][] scoringMatrix, int gapOpenPenalty, int gapExtendPenalty) {
        this.scoringMatrix = scoringMatrix;
        this.gapOpenPenalty = gapOpenPenalty;
        this.gapExtendPenalty = gapExtendPenalty;
    }

    @Override
    public int sequenceScore(UniqueSequence seq1, UniqueSequence seq2) throws DataException {
        return(fillDynamicMatrices(seq1, seq2).maxScore);
    }

    private FilledAlignmentMatrices fillDynamicMatrices(UniqueSequence seq1, UniqueSequence seq2) {
        int globalMaxScore = 0;
        int maxLine = 0;
        int maxColumn = 0;
        
        AlignmentMatrices matrices = initializeMatrices(seq1.getSequence().length + 1, seq2.getSequence().length + 1);
        int[][] scoreMatrix = matrices.getScoreMatrix();
        Direction[][] directionMatrix = matrices.getDirectionMatrix();
        
        for (int line = 1; line < seq1.getSequence().length + 1; line++){
            for(int column = 1; column < seq2.getSequence().length + 1; column++){
                
                int upPenalty;
                if (directionMatrix[line - 1][column] == Direction.UP) { //gap elonged
                    upPenalty = gapExtendPenalty;
                } else {
                    upPenalty = gapOpenPenalty;
                }
                
                int leftPenalty;
                if (directionMatrix[line][column - 1] == Direction.LEFT){ //gap elonged
                    leftPenalty = gapExtendPenalty;
                } else{
                    leftPenalty = gapOpenPenalty;
                }
                
                int upScore = scoreMatrix[line - 1][column] + upPenalty;
                int leftScore = scoreMatrix[line][column - 1] + leftPenalty;
                int diagonalScore = scoreMatrix[line - 1][column - 1] + scoringMatrix[seq1.getSequence()[line - 1]][seq2.getSequence()[column - 1]];
                
                int maxScore = Math.max(diagonalScore, (Math.max(upScore, leftScore)));
                
                if (maxScore < 0) {
                    scoreMatrix[line][column] = 0;
                    directionMatrix[line][column] = Direction.NOWHERE;
                } else{
                    scoreMatrix[line][column] = maxScore;
                    if (maxScore > globalMaxScore){
                        maxLine = line;
                        maxColumn = column;
                        globalMaxScore = maxScore;
                    }
                    if (maxScore == leftScore){
                        directionMatrix[line][column] = Direction.LEFT;
                    }
                    if (maxScore == upScore){
                        directionMatrix[line][column] = Direction.UP;
                    }
                    if (maxScore == diagonalScore){
                        directionMatrix[line][column] = Direction.DIAGONAL;
                    }
                }     
            }
        }
        return new FilledAlignmentMatrices(new AlignmentMatrices(scoreMatrix, directionMatrix), globalMaxScore, maxLine, maxColumn);
    }

    private static AlignmentMatrices initializeMatrices(int dim1, int dim2) {
        int[][] scoreMatrix = new int[dim1][dim2];
        Direction[][] directionMatrix = new Direction[dim1][dim2];

        for (int i = 1; i < dim1; i++) {
            scoreMatrix[i][0] = 0;
            directionMatrix[i][0] = Direction.UP;
        }
        for (int j = 1; j < dim2; j++) {
            scoreMatrix[0][j] = 0;
            directionMatrix[0][j] = Direction.LEFT;
        }
        return new AlignmentMatrices(scoreMatrix, directionMatrix);
    }

    private static class AlignmentMatrices {
        private final int[][] scoreMatrix;
        private final Direction[][] directionMatrix;

        public AlignmentMatrices(int[][] mainMatrix, Direction[][] directionMatrix) {
            this.scoreMatrix = mainMatrix;
            this.directionMatrix = directionMatrix;
        }

        public int[][] getScoreMatrix() {
            return scoreMatrix;
        }

        public Direction[][] getDirectionMatrix() {
            return directionMatrix;
        }
    }
    
    private static class FilledAlignmentMatrices{
        
        private final AlignmentMatrices matrices;
        private final int maxScore;
        private final int maxScoreLine;
        private final int maxScoreColumn;

        public FilledAlignmentMatrices(AlignmentMatrices matrices, int maxScore, int maxScoreLine, int maxScoreColumn) {
            this.matrices = matrices;
            this.maxScore = maxScore;
            this.maxScoreLine = maxScoreLine;
            this.maxScoreColumn = maxScoreColumn;
        }

        public AlignmentMatrices getMatrices() {
            return matrices;
        }

        public int getMaxScore() {
            return maxScore;
        }

        public int getMaxScoreLine() {
            return maxScoreLine;
        }

        public int getMaxScoreColumn() {
            return maxScoreColumn;
        }
        
    }
    
    

}
