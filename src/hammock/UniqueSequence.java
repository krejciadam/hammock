/*
 * Class represents single unique peptide sequence. Records all occurrencs of 
   this particualr sequence in dataset.
 */

package hammock;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Adam Krejci
 */
public class UniqueSequence implements Sizeable, Comparable<UniqueSequence>{
    
    private final int[] sequence;
    private final Map<String, Integer> labelsMap;
    private static final char[] aminoAcids = new char[]{'A', 'R', 'N', 'D', 'C',
                                                'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
                                                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 
                                                'V', 'B', 'Z', 'X', '*'};
    private static final Map<Character, Integer> nameToNum;
    
    static{
        Map<Character, Integer> NAMETONUM = new HashMap<>();
        for (int i = 0; i < aminoAcids.length; i++){
            NAMETONUM.put(aminoAcids[i], i);
        }
        nameToNum = Collections.unmodifiableMap(NAMETONUM);
    }
    
    /**
     * 
     * @param sequence Unique protein sequence (in standard AA letters, 
     *        upper or lower case)
     * @param labelsMap A map representing all sequence's labels.
     *          Key is label, value is count of occurrences of 
     *          this sequence with particular label
     */
    public UniqueSequence(String sequence, Map<String, Integer> labelsMap){
        this.labelsMap = labelsMap;
        this.sequence = new int[sequence.length()];
        String capitalSequence = sequence.toUpperCase();
        for (int i = 0; i < sequence.length(); i++){
            this.sequence[i] = nameToNum.get(capitalSequence.charAt(i));
        }
    }
    
    /**
     * Construction without labels map. Generates sequence with default label
     * "no_label" and count 1.
     * @param sequence Unique protein sequence (in standard AA letters, 
     *        upper or lower case)
     */
    public UniqueSequence (String sequence){
        Map<String, Integer> newLabelsMap = new HashMap<>();
        newLabelsMap.put("no_label", 1);
        this.labelsMap = newLabelsMap;
        this.sequence = new int[sequence.length()];
        String capitalSequence = sequence.toUpperCase();
        for (int i = 0; i < sequence.length(); i++){
            this.sequence[i] = nameToNum.get(capitalSequence.charAt(i));
        }
    }
    
    /**
     * Returns size of this sequence set, i.e. the number of occurrences 
     * of this particular sequence summed over all labels
     * @return Sum of occurrences of this sequence with all labels.
     */
    @Override
    public int size(){
        int size = 0;
        for (int partialSize : labelsMap.values()){
            size += partialSize;
        }
        return size;
    }

    /**
     * Returns number representation of instance's unique peptide string. 
     * @return unique paptide string as integer array.
     */
    public int[] getSequence() {
        return sequence;
    }
    
    /**
     * Returns String (standard AA letters) representation of instance's
     * unique peptide string.
     * @return Unique peptide string as String (standard AA letters, upper case)
     */
    public String getSequenceString() {
        char[] seq = new char[sequence.length];
        for (int i = 0; i < sequence.length; i++){
            seq[i] = aminoAcids[sequence[i]];
        }
        return new String(seq);
    }

    /**
     * Returns amino acid sequence. All operations using amino acids
     * (especially using a scoring matrix) assume amino acids in exactly 
     * this order. 
     * @return array of amino acids in correct order.
     */
    public static char[] getAminoAcidsAndSpecials() {
        return aminoAcids;
    }
    
    public static char[] getAminoAcids(){
        return Arrays.copyOfRange(aminoAcids, 0, 20);
    }
        
    /**
     * Returns a map representing occurrences of this sequence in the dataset. 
     *          Key is a label, value is count of occurrences of 
     *          this sequence with this label.
     * @return A map representing occurrences of sequences with particular labels.
     */
    public Map<String, Integer> getLabelsMap() {
        return labelsMap;
    }
    
    /**
     * Formates this object as a string to be written to a file.
     * @return 
     */
    public String getSavableString(){
        String res = this.getSequenceString();
        for (Map.Entry<String, Integer> entry : this.labelsMap.entrySet()){
            res = res.concat("|" + entry.getKey() + ":" + entry.getValue());
        }
        res = res.concat("\n");
        return res;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 79 * hash + Arrays.hashCode(this.sequence);
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
        final UniqueSequence other = (UniqueSequence) obj;
        return Arrays.equals(this.sequence, other.sequence);
    }

    /**
     * sequences are compared according to size, if size is the same, according to alphabetic order
     * @param o
     * @return 
     */
    @Override
        public int compareTo(UniqueSequence o) {
        if (o == this){
            return 0;
        }
        
        if (this.size() != o.size()){
            return this.size() - o.size();
        } else{
            return -(this.getSequenceString().compareTo(o.getSequenceString()));
        }
    }
        
        
  /*      ----------STATIC METHODS------------------*/
        
    public static List<UniqueSequence> sortSequences(List<UniqueSequence> sequences, String order) throws DataException{
        List<UniqueSequence> sortedList = sequences;
        switch(order){                
            case "size":
                Collections.sort(sortedList, Collections.reverseOrder(new UniqueSequenceSizeAlphabeticComparator()));
                break;
                
            case "alphabetic":
                Collections.sort(sortedList, Collections.reverseOrder(new UniqueSequenceAlphabeticComparator()));
                break;
                
            case "random":
                Collections.shuffle(sortedList, Hammock.random);
                break;
            
            default:
                throw new DataException("Incorrect sequence order defined. Use one of: size, alphabetic, random");
        }
        return sortedList;
    }    
        
        
}

/**
 * Compares sequences according to sum of counts of all labels. In case of
 * equality, compares alphabetically.
 *
 * @author Adam Krejci
 */
class UniqueSequenceSizeAlphabeticComparator implements Comparator<UniqueSequence> {

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        int result = new SizeComparator().compare(o1, o2);
        if (result == 0) {
            result = new UniqueSequenceAlphabeticComparator().compare(o1, o2);
        }
        return result;
    }
}

/**
 * Compares sequences alphabetically
 *
 * @author Adam Krejci
 */
class UniqueSequenceAlphabeticComparator implements Comparator<UniqueSequence> {

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        return o1.getSequenceString().compareTo(o2.getSequenceString());
    }
}
