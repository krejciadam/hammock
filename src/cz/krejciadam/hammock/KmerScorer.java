/*
 * The class scores sequences on the basis of kmer ocurrences
 */
package cz.krejciadam.hammock;
 
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
 
/**
 *
 * @author Adam Krejci
 */
public class KmerScorer implements SequenceScorer {
 
    private final int k;
    private final Alphabet alphabet;
 
    public KmerScorer(int k, Alphabet alphabet) {
        this.k = k;
        this.alphabet = alphabet;
    }
    
    @Override
    public int sequenceScore(UniqueSequence seq1, UniqueSequence seq2) throws DataException {
        Set<List<Integer>> kmers1 = getAllKmers(seq1);
        Set<List<Integer>> kmers2 = getAllKmers(seq2);
        kmers1.retainAll(kmers2);
        return(kmers1.size());
    }
 
 
    public Set<List<Integer>> getAllKmers(UniqueSequence seq) {
        Set<List<Integer>> res = new HashSet<>();
        Integer[] fullSeq = alphabet.toSequence(seq.getSequenceString());
        for (int i = 0; i < (seq.getSequenceString().length() - k + 1); i++) {
            List<Integer> kmer = Collections.unmodifiableList(Arrays.asList(Arrays.copyOfRange(fullSeq, i, i + k)));
            res.add(kmer);
        }
        return (res);
    }
}