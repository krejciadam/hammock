/*
 * Class represents cluster of sequences. Cluster contains any number of unique
 * peptide sequences. Cluster has its own unique id. Users should guarantee
 * id's uniqueness. Cluster may have MSA, HMM and HH (= hmm for HMM-HMM comparison)
 constructed. 
 */
package hammock;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Adam Krejci
 */
public class Cluster implements Sizeable, Comparable<Cluster> {

    private final List<UniqueSequence> sequences;
    private final int id;
    private int size;
    private final Map<String, Integer> labelsMap;
    private boolean hasMSA = false;
    private boolean hasHMM = false;
    private boolean hasHH = false;

    public Cluster(Collection<UniqueSequence> sequences, int id) {
        this.sequences = new ArrayList<>(sequences);
        this.id = id;
        int sizeSum = 0;
        labelsMap = new HashMap<>();
        for (UniqueSequence seq : sequences) {
            sizeSum += seq.size();
            updateLabelsMap(seq);
        }
        this.size = sizeSum;
    }

    /**
     * Inserts another UniqueSequence into the cluster. Immediately after this
     * operation, cluster's hasMSA(), hasHMM() and hasHH() will return false.
     *
     * @param sequence A sequence to be added to the cluster.
     * @throws hammock.DataException
     */
    public void insert(UniqueSequence sequence) throws DataException {
        if (sequences.contains(sequence)) {
            throw new DataException("Trying to insert unique sequence "
                    + sequence.getSequenceString() + " into cluster "
                    + this.id + ", which already contains this sequence. ");
        }
        sequences.add(sequence);
        this.size += sequence.size();
        updateLabelsMap(sequence);

        this.hasMSA = false;
        this.hasHMM = false;
        this.hasHH = false;
    }
    
    /**
     * Iserts multiple sequences into the cluster.
     * @param sequences sequences to be inserted
     * @throws DataException 
     */
    public void insertAll(Collection<UniqueSequence> sequences) throws DataException{
        for (UniqueSequence seq : sequences){
            insert(seq);
        }
    }

    private void updateLabelsMap(UniqueSequence seq) {
        for (Map.Entry<String, Integer> entry : seq.getLabelsMap().entrySet()) {
            Integer count = labelsMap.get(entry.getKey());
            if (count == null){
                count = 0;
            }
            count += entry.getValue();
            labelsMap.put(entry.getKey(), count);
        }
    }

    /**
     * Returns the vector of label counts. Uses the global labels list from the main class
     * @return 
     */
    public int[] getLabelCountVector() {
        int[] vector = new int[Hammock.getLabels().size()];
        int i = 0;
        for (String label : Hammock.getLabels()) {
            Integer count = labelsMap.get(label);
            if (!(count == null)){
                vector[i] = labelsMap.get(label);
            } else{
                vector[i] = 0;
            }
            i++;
        }
        return (vector);
    }

    /**
     * Returns the number of unique sequences in this cluster, i.e. every sequence
     * present counted as 1, no matter what are the occurrences of this sequence
     * with different labels.
     *
     * @return Number of unique sequence strings in this cluster.
     */
    public int getUniqueSize() {
        return sequences.size();
    }

    public List<UniqueSequence> getSequences() {
        return sequences;
    }

    public int getId() {
        return id;
    }

    public boolean hasMSA() {
        return hasMSA;
    }

    public boolean hasHMM() {
        return hasHMM;
    }

    public boolean hasHH() {
        return hasHH;
    }

    public void setAsHasMSA() throws IOException {
        this.hasMSA = true;
    }

    public void setAsHasHMM() {
        this.hasHMM = true;
    }

    public void setAsHasHH() {
        this.hasHH = true;
    }

    /**
     * Returns the sum of sizes of UniqueSequence instances contained in this
     * cluster, i.e. the size with respect to occurrences with labels.
     *
     * @return
     */
    @Override
    public int size() {
        return size;
    }


    /**
     * Returns cluster's unique sequences in form of amino acid sequences in
     * fasta format. Sequences are returned in reversed natural order.
     *
     * @return
     */
    public String getFastaString() {
        int seqId = 1;
        String res = "";
        Collections.sort(sequences, Collections.reverseOrder());
        for (UniqueSequence sequence : sequences) {
            res = res.concat(">" + id + "_" + seqId + "\n" + sequence.getSequenceString() + "\n");
            seqId++;
        }
        return res;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 79 * hash + this.id;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Cluster other = (Cluster) obj;
        return this.id == other.id;
    }

    @Override
    public int compareTo(Cluster o) {
        if (size != o.size()) {
            return size - o.size();
        } else {
            return id - o.getId();
        }
    }

}
