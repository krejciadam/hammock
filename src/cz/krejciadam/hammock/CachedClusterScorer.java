/*
 * Performs fast cached cluser scoring. Currently only supports complete-linkage.
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Adam Krejci
 */

/**
 * Warning! It is needed to call the join() method appropriately when this scorer is used.
 * @author Adam Krejci
 */
public class CachedClusterScorer implements ClusterScorer {
    
    private final Map<Integer, Integer> indexMap; //Index. The id of each cluster in the matrix gets a matrix x-coordinate.
    private final DynamicMatrix matrix;
    private final int sizeLimit;
    private final ClusterScorer defaultClusterScorer;
    
    public CachedClusterScorer(ClusterScorer defaultClusterScorer, int sizeLimit){
        this.sizeLimit = sizeLimit;
        this.defaultClusterScorer = defaultClusterScorer;
        indexMap = Collections.synchronizedMap(new HashMap<Integer, Integer>());
        matrix = new DynamicMatrix();
    }
    
    @Override
    public int clusterScore(Cluster cl1, Cluster cl2) throws DataException{
        if (cl1.getUniqueSize() < sizeLimit || cl2.getUniqueSize() < sizeLimit){
            return(defaultClusterScorer.clusterScore(cl1, cl2));
        }
        Integer score;
        if (indexMap.containsKey(cl1.getId())) {    //both clusters are in the index
            if (indexMap.containsKey(cl2.getId())) {
                Integer index1 = indexMap.get(cl1.getId());
                Integer index2 = indexMap.get(cl2.getId());
                score = matrix.get(index1, index2);
                if (score != null) {
                    return score;
                } else {
                    score = defaultClusterScorer.clusterScore(cl1, cl2);
                    matrix.set(index1, index2, score);
                }
            } else { //cl1 is and cl2 is not
                int matrixIndex = matrix.getMatrixIndex(); //prepare new row
                score = defaultClusterScorer.clusterScore(cl1, cl2);
                int firstIndex = indexMap.get(cl1.getId());
                matrix.add(firstIndex, matrixIndex, score);
                indexMap.put(cl2.getId(), matrixIndex);
            }
        } else {
            if (indexMap.containsKey(cl2.getId())) { //cl2 is and cl1 is not
                int matrixIndex = matrix.getMatrixIndex(); //prepare new row
                score = defaultClusterScorer.clusterScore(cl1, cl2);
                int firstIndex = indexMap.get(cl2.getId());
                matrix.add(firstIndex, matrixIndex, score);
                indexMap.put(cl1.getId(), matrixIndex);
            } else {                                //neither is
                int firstIndex = matrix.getMatrixIndex(); //prepare new row
                int secondIndex = matrix.getMatrixIndex();
                score = defaultClusterScorer.clusterScore(cl1, cl2);
                matrix.addEmpty(firstIndex);
                matrix.addEmpty(secondIndex);
                indexMap.put(cl1.getId(), firstIndex);
                indexMap.put(cl2.getId(), secondIndex);
            }
        }
        return score;
    }
    
    
    public synchronized void join(Cluster cl1, Cluster cl2, int newId){
        /*Clusters too small for saving - nothing to merge*/
        if (cl1.getUniqueSize() < sizeLimit || cl2.getUniqueSize() < sizeLimit){
            return;
        }
        
        /*An attempt to merge two non-ready rows cannont happen, if "join" is not called after all 
        threads using "add" finished*/
        
        if (indexMap.containsKey(cl1.getId())) {
            int index1 = indexMap.get(cl1.getId());
            if (indexMap.containsKey(cl2.getId())) {
                int index2 = indexMap.get(cl2.getId());
                List<Integer> line1 = matrix.get(index1);
                List<Integer> line2 = matrix.get(index2);
                List<Integer> newLine = new ArrayList<>();
                for (int i = 0; i < line1.size(); i++) {
                    Integer num1 = line1.get(i);
                    Integer num2 = line2.get(i);
                    if (num1 != null && num2 != null) {
                        newLine.add(i, Math.min(line1.get(i), line2.get(i))); //calculate the score of the merged cluster to all others
                    } else {
                        newLine.add(i, null);
                    }
                }
                matrix.remove(index1); //remove old rows
                matrix.remove(index2);
                int newIndex = matrix.getMatrixIndex(); //pripare for row insertion
                matrix.add(newLine, newIndex);
                indexMap.remove(cl1.getId()); //remove old indexes
                indexMap.remove(cl2.getId());
                indexMap.put(newId, newIndex);
            } else {    //no values for both clusters saved
                matrix.remove(index1);
                indexMap.remove(cl1.getId());
            }
        } else { //no values for both clusters saved
            if (indexMap.containsKey(cl2.getId())) {
                int index2 = indexMap.get(cl2.getId());
                matrix.remove(index2);
                indexMap.remove(cl2.getId());
            }
        }        
    }    
}

class DynamicMatrix {

    private final List<List<Integer>> matrix;
    private final Set<Integer> dirty;

    public DynamicMatrix() {
        matrix = Collections.synchronizedList(new ArrayList<List<Integer>>());
        dirty = Collections.synchronizedSet(new HashSet<Integer>());

        /* Fill first two rows (and "columns") of the matrix. Row with index of 
        2 is ready for new valuest*/
        List<Integer> list = Collections.synchronizedList(new ArrayList<Integer>());
        list.add(null);
        matrix.add(list);
        list = Collections.synchronizedList(new ArrayList<Integer>());
        list.add(null);
        list.add(null);
        matrix.add(list);
        dirty.add(1);
    }

    /**
     * Returns row (=column) with the index of lineIndex
     * @param rowIndex coordinates of the row (=columns)
     * @return row saved on lineIndex
     */
    public List<Integer> get(int rowIndex) {
        List<Integer> line = new ArrayList<>(matrix.get(rowIndex));
        for (int i = rowIndex + 1; i < matrix.size(); i++) {
            line.add(matrix.get(i).get(rowIndex));
        }
        return new ArrayList<>(line);
    }

    /**
     * Returns matrix value on coordinates
     * @param rowIndex x-coordinate
     * @param columnIndex y-coordinate
     * @return value placed in the matrix on index1,index2 coordinates
     */
    public Integer get(int rowIndex, int columnIndex) {
        List<Integer> line = matrix.get(Math.max(rowIndex, columnIndex));
        if (line != null) {
            return line.get(Math.min(rowIndex, columnIndex));
        } else {
            return null;
        }
    }

    /**
     * Sets matrix value on specific coordinates
     *
     * @param rowIndex x-coordinate
     * @param columnIndex y-coordinate
     * @param value value to be set
     */
    public void set(int rowIndex, int columnIndex, int value) {
        int min = Math.min(rowIndex, columnIndex);
        int max = Math.max(rowIndex, columnIndex);
        List<Integer> line = matrix.get(max);
        if (line != null) {
            line.set(min, value);
            matrix.set(max, line);
        } else {
            //Do nothing. The entire row is null, which means another thread is 
            //generating it right now. Writing here would mean the other thread
            //will override the value.
        }
    }

    /**
     * Marks a row (=column) in the matrix as no longer used.
     * @param rowIndex x-coordinate
     */
    public void remove(int rowIndex) {
        dirty.add(rowIndex);
    }

    /**
     * Inserts an empty row into the matrix
     *
     * @param matrixIndex the row coordinate to put the empty line on. It is NECESSARY
     * that this value was generated by getMatrixIndex().
     * @throws PhageException
     */
    public void addEmpty(int matrixIndex){
        List<Integer> line = new ArrayList<>();
        for (int i = 0; i < matrixIndex + 1; i++) {
            line.add(null);
        }
        this.add(line, matrixIndex);
    }

    /**
     * Inserts a row containing one value, otherwise empty.
     *
     * @param index index of the value in the row
     * @param rowIndex x-coordinate of the row in the matrix. It is NECESSARY that
     * this value was generated using getmatrixIndex();
     * @param value The value to save.
     * @throws PhageException
     */
    public void add(int index, int rowIndex, int value){
        List<Integer> line = new ArrayList<>();
        for (int i = 0; i < matrix.size(); i++) {
            line.add(null);
        }
        line.set(index,value);
        this.add(line, rowIndex);
    }

    /**
     * Inserts a row on a specific position in the matrix
     *
     * @param row The row to insert
     * @param rowIndex x-coordinate of the row in the matrix. It is NECESSARY that
     * this value was generated using getmatrixIndex();
     * @throws PhageException
     */
    public void add(List<Integer> row, int rowIndex){
        List<Integer> firstPart = row.subList(0, rowIndex + 1);
        matrix.set(rowIndex, firstPart);
        for (int i = rowIndex + 1; i < matrix.size(); i++) {
            List<Integer> currentLine = matrix.get(i);
            if ((currentLine != null) && (i < row.size())) {
                currentLine.set(rowIndex, row.get(i));
                matrix.set(i, currentLine);
            } else {
                //Do nothing. We would write null, but it is already here or soon will be.
            }

        }
    }

    /**
     * Returns matrix x-coordinate where a new row can be inserted. It is necessary
     * to call this method before every add() or addEmpty()
     *
     * @return matrix x-coordinate to put a new row at
     */
    public synchronized int getMatrixIndex() {
        int matrixIndex;
        if (!dirty.isEmpty()) {
            matrixIndex = dirty.iterator().next();
            dirty.remove(matrixIndex);
        } else {
            matrixIndex = matrix.size();
            matrix.add(null);
        }
        return matrixIndex;
    }
    
}
