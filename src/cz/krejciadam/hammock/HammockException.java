/*
 * Project's own exception
 */

package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class HammockException extends Exception {
    
    public HammockException(String s){
        super(s);
    }
    
    public HammockException(Throwable cause){
        super(cause);
    }
    
    public HammockException(String s, Throwable cause){
        super(s, cause);
    }
    
}
