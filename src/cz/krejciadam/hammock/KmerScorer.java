/*
 * The class serves to score kmers on the basis of their common ocurrence in sequences
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Adam Krejci
 */
public class KmerScorer implements SequenceScorer {

    private final int k;
    private final ScoreMap scoreMap;
    private final Alphabet alphabet;

    public KmerScorer(int k, Collection<UniqueSequence> sequences, Alphabet alphabet) {
        this.k = k;
        this.scoreMap = new ScoreMap(k, sequences, alphabet);
        this.alphabet = alphabet;
    }

    @Override
    public int sequenceScore(UniqueSequence seq1, UniqueSequence seq2) throws DataException {
        if (seq1.getSequence().length != k || seq2.getSequence().length != k) {
            throw new DataException("Error. Trying to use a kmer scorer intitialized for k = " + k + " on"
                    + "sequences of lnegths " + seq1.getSequence().length + " and " + seq2.getSequence().length + " .");
        }
        int index1 = kmer2index(alphabet.toSequence(seq1.getSequenceString()), alphabet.getAlphabetSize());
        int index2 = kmer2index(alphabet.toSequence(seq2.getSequenceString()), alphabet.getAlphabetSize());
        return(scoreMap.get(index1, index2));
    }

    private int kmer2index(int[] kmer, int alphabetLength) {
        int res = 0;
        for (int i = 0; i < kmer.length; i++) {
            res += kmer[i] * (int) Math.pow(alphabetLength, (kmer.length - i - 1));
        }
        return (res);
    }

    private List<int[]> getAllKmers(UniqueSequence seq, int k, Alphabet alphabet) {
        List<int[]> res = new ArrayList<>();
        int[] fullSeq = alphabet.toSequence(seq.getSequenceString());
        for (int i = 0; i < (seq.getSequenceString().length() - k + 1); i++) {
            int[] kmer = Arrays.copyOfRange(fullSeq, i, i + k);
            res.add(kmer);
        }
        return (res);
    }

    class ScoreMap {

        private final int[][] map;

        public ScoreMap(int k, Collection<UniqueSequence> sequences, Alphabet alphabet) {
            map = initializeScoreMap(k, alphabet.getAlphabetSize());
            fillScoreMap(k, sequences, alphabet);
        }

        private int[][] initializeScoreMap(int k, int alphabetLength) {
            int[][] scoreMap;
            int kmerCount = (int) Math.pow(alphabetLength, k);
            scoreMap = new int[kmerCount][];
            for (int i = 0; i < kmerCount; i++) {
                scoreMap[i] = new int[kmerCount - i];
                for (int j = 0; j < (kmerCount - i); j++) {
                    scoreMap[i][j] = 0;
                }
            }
            return (scoreMap);
        }

        private void fillScoreMap(int k, Collection<UniqueSequence> sequences, Alphabet alphabet) {
            for (UniqueSequence seq : sequences) {
                Set<Integer> indicesSet = new HashSet<>(); //ensure we only count unique kmers
                for (int[] kmer : getAllKmers(seq, k, alphabet)) {
                    indicesSet.add(kmer2index(kmer, alphabet.getAlphabetSize()));
                }
                List<Integer> indicesList = new ArrayList<>();
                indicesList.addAll(indicesSet);
                for (int i = 0; i < (indicesList.size()); i++) {
                    for (int j = i; j < indicesList.size(); j++) { //not j + 1. i == j means a single kmer and we want to count single kmer ocurrences
                        int[] indices = getCorrectIndices(indicesList.get(i), indicesList.get(j));
                        map[indices[0]][indices[1]] += 1;
                    }
                }
            }
        }

        public int get(int index1, int index2) {
            int[] indices = getCorrectIndices(index1, index2);
            return (map[indices[0]][indices[1]]);
        }


        private int[] getCorrectIndices(int index1, int index2) {
            if (index1 > index2) {
                int pom = index1;
                index1 = index2;
                index2 = pom;
            }
            return (new int[]{index1, index2 - index1});
        }
    }

}
