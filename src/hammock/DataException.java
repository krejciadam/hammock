/*
 * This exception reports errors in data structures
 */

package hammock;

/**
 *
 * @author Adam Krejci
 */
public class DataException extends HammockException{

    public DataException(String s) {
        super(s);
    }

    public DataException(Throwable cause) {
        super(cause);
    }

    public DataException(String s, Throwable cause) {
        super(s, cause);
    }
    
}
