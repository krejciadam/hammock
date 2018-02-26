/* 
 * An amino acid alphabet
 */
package cz.krejciadam.hammock;
 
import java.util.HashMap;
import java.util.Map;
 
/**
 *
 * @author Adam Krejci
 */
public abstract class Alphabet {
    
    Map<Character, Integer> aaMap;
    int alphabetSize;
    
    public Integer[] toSequence(String sequenceString){
        Integer[] res = new Integer[sequenceString.length()];
        for (int i = 0; i < sequenceString.length(); i++){
            res[i] = aaMap.get(sequenceString.charAt(i));
        }
        return(res);
    }
    
    
    public static Map<Character, Integer> getAAMap(String[] groups){
        Map<Character, Integer> res = new HashMap<>();
        for (int i = 0; i < groups.length; i++){
            String s = groups[i];
            for (int j = 0; j < s.length(); j++){
                res.put(s.charAt(j), i);
            }
        }
        return(res);
    }
    
    
     /**
     * Returns an array containing all the alphabet letters as integers.
     * @return 
     */
    public int getAlphabetSize(){
        return(alphabetSize);
    }
    
    
}