/**
 * Represents a pair of Comparable objects. Order of these objects does
 * not matter, bigger object is always first, no matter what order were 
 * objects passed to constructor.
 * @author Adam Krejci
 * @param <T> Type, must extend Comparable. Both members must be of the same type
 */

package cz.krejciadam.hammock;

import java.util.Objects;

/**
 *
 * @author Adam Krejci
 * @param <T>
 */
public class UnorderedPair<T extends Comparable<? super T>>{
    private final T smaller;
    private final T bigger;

    public UnorderedPair(T first, T second) {
        if (first.compareTo(second) >= 0){
            bigger = first;
            smaller = second;
        } else{
            bigger = second;
            smaller = first;
        }
    }

    public T getSmaller() {
        return smaller;
    }

    public T getBigger() {
        return bigger;
    }  

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 83 * hash + Objects.hashCode(this.smaller);
        hash = 83 * hash + Objects.hashCode(this.bigger);
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
        final UnorderedPair<?> other = (UnorderedPair<?>) obj;
        if (!Objects.equals(this.smaller, other.smaller)) {
            return false;
        }
        return Objects.equals(this.bigger, other.bigger);
    }
}