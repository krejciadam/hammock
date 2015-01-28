/*
 * This exception reports reading of incorrectly formatted file
 */

package hammock;

/**
 *
 * @author Adam Krejci
 */
public class FileFormatException extends HammockException{

    public FileFormatException(String s) {
        super(s);
    }

    public FileFormatException(Throwable cause) {
        super(cause);
    }

    public FileFormatException(String s, Throwable cause) {
        super(s, cause);
    }
    
}
