/*
 * Compares Sizeable objects according to their size
 */

package hammock;

import java.util.Comparator;

/**
 *
 * @author Adam Krejci
 */
public class SizeComparator implements Comparator<Sizeable>{
    @Override
    public int compare(Sizeable o1, Sizeable o2) {
        int size1 = o1.size();
        int size2 = o2.size();

        return size1 -size2;
    } 
}
