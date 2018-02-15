/*
 * Reduced amino acid alphabet SDM12
 * Singleton. 
 */
package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class AlphabetSdm12 extends Alphabet {
    
    private static final AlphabetSdm12 INSTANCE = new AlphabetSdm12();

    private AlphabetSdm12() {
        super.aaMap = Alphabet.getAAMap(new String[]{
                                   "A", "D", "KER", 
                                   "N", "TSQ", "YF",
                                   "LIVM", "C", "W", 
                                   "H", "G", "P"});
        super.alphabetSize = 12;
    }

    public static AlphabetSdm12 getInstance(){
        return INSTANCE;
    }
}
