/*
 * Performs fast cached cluser scoring. Currently only supports complete-linkage.
 */
package hammock;

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
 * @author akrejci
 */
public class CachedClusterScorer implements ClusterScorer {
    
    private static Map<Integer, Integer> indexMap; //Index. Kazdemu id clusteru v matici prirazuje index odpovidajiciho radku matice.
    private static DynamicMatrix matrix;
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
        if (indexMap.containsKey(cl1.getId())) {    //oba clustery jsou v indexu
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
            } else { //cl1 je a cl2 neni
                int matrixIndex = matrix.getMatrixIndex(); //priprava noveho radku
                score = defaultClusterScorer.clusterScore(cl1, cl2);
                int firstIndex = indexMap.get(cl1.getId());
                matrix.add(firstIndex, matrixIndex, score);
                indexMap.put(cl2.getId(), matrixIndex);
            }
        } else {
            if (indexMap.containsKey(cl2.getId())) { //cl2 je a cl1 neni
                int matrixIndex = matrix.getMatrixIndex(); //priprava noveho radku
                score = defaultClusterScorer.clusterScore(cl1, cl2);
                int firstIndex = indexMap.get(cl2.getId());
                matrix.add(firstIndex, matrixIndex, score);
                indexMap.put(cl1.getId(), matrixIndex);
            } else {                                //neni  ani jeden
                int firstIndex = matrix.getMatrixIndex(); //priprava noveho radku
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
        /*Prilis male clustery pro ukladani - nemame co slucovat*/
        if (cl1.getUniqueSize() < sizeLimit || cl2.getUniqueSize() < sizeLimit){
            return;
        }
        
        /*Pokus o slouceni dvou jeste nehotovych radku nemuze nastat, pokud join nevolame po ukonceni vsech vlaken pouzivajicich add*/
        
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
                        newLine.add(i, Math.min(line1.get(i), line2.get(i))); //vypocet skore slouceneho clusteru ke vsem ostatnim
                    } else {
                        newLine.add(i, null);
                    }
                }
                matrix.remove(index1); //odstranime stare radky
                matrix.remove(index2);
                int newIndex = matrix.getMatrixIndex(); //priprava vlozeni radku
                matrix.add(newLine, newIndex);
                indexMap.remove(cl1.getId()); //odstraneni starych indexu
                indexMap.remove(cl2.getId());
                indexMap.put(newId, newIndex);
            } else {    //nemame ulozeny hodnoty pro oba clustery
                matrix.remove(index1);
                indexMap.remove(cl1.getId());
            }
        } else { //nemame ulozeny hodnoty pro oba clustery
            if (indexMap.containsKey(cl2.getId())) {
                int index2 = indexMap.get(cl2.getId());
                matrix.remove(index2);
                indexMap.remove(cl2.getId());
            }
        }        
    }    
}

class DynamicMatrix {

    private static List<List<Integer>> matrix;
    private static Set<Integer> dirty;

    public DynamicMatrix() {
        matrix = Collections.synchronizedList(new ArrayList<List<Integer>>());
        dirty = Collections.synchronizedSet(new HashSet<Integer>());

        /* Naplneni prvnich dvou radku (a "sloupcu") matice. Radek s indexem 2
         je pripraven na zapis hodnot*/
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
     * Vrati radek (= sloupec) se zadanym indexem.
     *
     * @param lineIndex souradnice radku (= sloupce)
     * @return radek, ktery je v matici ulozen na zadanem indexu
     */
    public List<Integer> get(int lineIndex) {
        List<Integer> line = new ArrayList<>(matrix.get(lineIndex));
        for (int i = lineIndex + 1; i < matrix.size(); i++) {
            line.add(matrix.get(i).get(lineIndex));
        }
        return new ArrayList<>(line);
    }

    /**
     * Vrati cislo z matice na zadanych souradnicich
     *
     * @param index1 prvni souradnice
     * @param index2 druha souradnice
     * @return cislo, ktere je v matici na zadanych souradnicich
     */
    public Integer get(int index1, int index2) {
        List<Integer> line = matrix.get(Math.max(index1, index2));
        if (line != null) {
            return line.get(Math.min(index1, index2));
        } else {
            return null;
        }
    }

    /**
     * Nastavi hodnotu v matici na zadanych souradnicich
     *
     * @param lineIndex prvni souradnice
     * @param columnIndex druha souradnice
     * @param value nastavovana hodnota
     */
    public void set(int lineIndex, int columnIndex, int value) {
        int min = Math.min(lineIndex, columnIndex);
        int max = Math.max(lineIndex, columnIndex);
        List<Integer> line = matrix.get(max);
        if (line != null) {
            line.set(min, value);
            matrix.set(max, line);
        } else {
            //nic neudelame. cela radka je null, tedy jine vlakno ji zrovna generuje. Kdybychom 
            //ted radku pripsali, jine vlakno nam ji hned prepise.
        }
    }

    /**
     * Odstrani radek (=sloupec) z matice, resp. oznaci jej z dale nepouzivany
     *
     * @param lineIndex souradnice radku k odstraneni
     */
    public void remove(int lineIndex) {
        dirty.add(lineIndex);
    }

    /**
     * Vlozi do matice prazdny radek
     *
     * @param matrixIndex souradnice, kam ma byt radek vlozen. Je NUTNE, aby
     * tato hodnota byla predem ziskana pomoci getMatrixIndex();
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
     * Vlozi do matice radek s jednou hodnotou, jinak prazdny.
     *
     * @param index souradnice hodnoty v radku
     * @param matrixIndex souradnice radku v matici. Je NUTNE, aby tato hodnota
     * byla predem ziskana pomoci getMatrixIndex();
     * @param value ukladana hodnota.
     * @throws PhageException
     */
    public void add(int index, int matrixIndex, int value){
        List<Integer> line = new ArrayList<>();
        for (int i = 0; i < matrix.size(); i++) {
            line.add(null);
        }
        line.set(index,value);
        this.add(line, matrixIndex);
    }

    /**
     * Vlozi radek na zadanou pozici v matici.
     *
     * @param line Radek k vlozeni
     * @param index souradnice radku v matici. Je NUTNE, aby tato hodnota byla
     * predem ziskana pomoci getMatrixIndex();
     * @throws PhageException
     */
    public void add(List<Integer> line, int index){
        List<Integer> firstPart = line.subList(0, index + 1);
        matrix.set(index, firstPart);
        for (int i = index + 1; i < matrix.size(); i++) {
            List<Integer> currentLine = matrix.get(i);
            if ((currentLine != null) && (i < line.size())) {
                currentLine.set(index, line.get(i));
                matrix.set(i, currentLine);
            } else {
                //nebudeme delat nic. V tomto pripade bychom zapsali null, ale ten uz na miste bud je, nebo brzy bude
            }

        }
    }

    /**
     * Vrati souradnici radku v matici, kam ma byt vlozen dalsi vkladany radek.
     * Tuto metodu je nutno volat pred kazdym vkladanim do matice, tedy kazdym
     * volanim add() nebo addEmpty()
     *
     * @return souradnice v matici kam muze byt vlozen dalsi radek
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
