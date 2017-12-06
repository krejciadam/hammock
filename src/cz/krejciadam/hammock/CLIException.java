/*
 * An exception for command line argument errors.
 */
package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class CLIException extends HammockException{

    public CLIException(String s) {
        super(s);
    }

    public CLIException(Throwable cause) {
        super(cause);
    }

    public CLIException(String s, Throwable cause) {
        super(s, cause);
    }
    
}