/*
 * Class performs greedy clustering. 
 */

package hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class GreedySequenceClusterer implements SequenceClusterer {
    
    private final int threshold;

    
    /**
     * Constructor.
     * @param threshold Minimum score for a sequence to be added into cluster
     *                  during greedy clustering.
     */
    public GreedySequenceClusterer(int threshold) {
        this.threshold = threshold;
    }

    /**
     * Performs greedy sequence clustering. First, sequences are ordered. 
     * Then fist sequence is chosen as pivot and all sequences similar enough
     * (above threshold) to pivot are selected. Pivot and selected sequences
     * are reported as a cluster and removed from list. First sequence from
     * list is chosen as pivot and the cycle repeats.
     * 
     * This method works in O(n^2) and 
     * is deterministic, though determinism is achieved only by (inner) ordering
     * of input sequences by size and alphabet. 
     * @param sequences
     * @param scorer
     * @return
     * @throws InterruptedException
     * @throws ExecutionException 
     */
    @Override
    public List<Cluster> cluster(Collection<UniqueSequence> sequences, ShiftedScorer scorer) throws InterruptedException, ExecutionException{
        int nThreads = Hammock.nThreads;
        List<Cluster> result = new ArrayList<>();
        List<UniqueSequence> sortedList = new ArrayList<>(sequences);
        Collections.sort(sortedList, Collections.reverseOrder(new UniqueSequenceSizeAlphabeticComparator()));
        int id = 1;
        CompletionService<ListSearchResult> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        List<ListSearchResult> resultList;
        while (sortedList.size() > 0){
            UniqueSequence pivot = sortedList.get(0);
            sortedList.remove(0);
            int forOneThread = sortedList.size() / nThreads; //while number division
            int i = 1;
            for (; i < nThreads; i++){
                List<UniqueSequence> subList = sortedList.subList((i - 1) * forOneThread, i * forOneThread);
                resultPool.submit(new ListSearchRunner(subList, pivot, threshold, scorer, i));
            }
            List<UniqueSequence> subList = sortedList.subList((i - 1) * forOneThread, sortedList.size()); //last piece of list
            resultPool.submit(new ListSearchRunner(subList, pivot, threshold, scorer, i));  
            resultList = new ArrayList<>();
            for (int j = 0; j < nThreads; j++){
                resultList.add(resultPool.take().get());                
            }
            Collections.sort(resultList);
            sortedList = new ArrayList<>();
            
            List<UniqueSequence> clusterSequences = new ArrayList<>();
            clusterSequences.add(pivot);
            for (ListSearchResult res : resultList){
                sortedList.addAll(res.getOtherSequences());
                for (UniqueSequence seq : res.getSelectedSequences()){
                    clusterSequences.add(seq);
                }
            }     
            Cluster cl = new Cluster(clusterSequences, id);
            result.add(cl); 
            id ++;
            if (nThreads >= sortedList.size() + 2){
                nThreads = 1;
            } 
        }
        return result;
    }
}

/**
 * Helper class for parallel processing. When searching through a list 
 * of sequences with a pivot sequences, we divide the list in two parts - 
 * sequences similar enough to pivot (selected) and the other sequences.
 * @author Adam Krejci
 */
class ListSearchResult implements Comparable<ListSearchResult>{
    private final List<UniqueSequence> selectedSequences;
    private final List<UniqueSequence> otherSequences;
    private final int order;

    public ListSearchResult(List<UniqueSequence> selectedSequences, List<UniqueSequence> otherSequences, int order) {
        this.selectedSequences = selectedSequences;
        this.otherSequences = otherSequences;
        this.order = order; //we remember the order, so we do not need to sort.
    }

    public List<UniqueSequence> getSelectedSequences() {
        return selectedSequences;
    }

    public List<UniqueSequence> getOtherSequences() {
        return otherSequences;
    }

    public int getOrder() {
        return order;
    }
    @Override
    public int compareTo(ListSearchResult o) {
        return this.getOrder() - o.getOrder();
    }
}


/**
 * Helper class for parallel processing. Searches through a given list of sequences
 * using a pivot sequences. All sequences similar to pivot above given threshild
 * are reported as selected, all the other sequences are reported as other sequences.
 * 
 * @author Adam Krejci
 */
class ListSearchRunner implements Callable<ListSearchResult>{
    
    private final List<UniqueSequence> sequenceList;
    private final UniqueSequence pivot;
    private final int scoreThreshold;
    private final Scorer scorer;
    private final int order; //pamatujeme si poradi, abychom zabranili nutnosti trideni

    public ListSearchRunner(List<UniqueSequence> phageList, UniqueSequence pivot, int scoreThreshold, Scorer scorer, int order) {
        this.sequenceList = phageList;
        this.pivot = pivot;
        this.scoreThreshold = scoreThreshold;
        this.scorer = scorer;
        this.order = order;
    }
    
    @Override
    public ListSearchResult call() throws DataException{
        int score;
        List<UniqueSequence> selectedSequences = new ArrayList<>();
        List<UniqueSequence> otherSequences = new ArrayList<>();
        for (UniqueSequence comparedSequence : sequenceList){
            score = scorer.score(pivot, comparedSequence);
            if (score >= scoreThreshold){
                selectedSequences.add(comparedSequence);
            } else{
                otherSequences.add(comparedSequence);
            }
        }
        return new ListSearchResult(selectedSequences, otherSequences, order);
    }
}

/**
 * Compares sequences according to sum of counts of all labels. In case 
 * of equality, compares alphabetically.
 * @author Adam Krejci
 */
class UniqueSequenceSizeAlphabeticComparator implements Comparator<UniqueSequence>{

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        int result = new SizeComparator().compare(o1, o2);
        if (result == 0){
            result = new UniqueSequenceAlphabeticComparator().compare(o1, o2);
        }
        return result;
    }  
}

/**
 * Compares sequences alphabetically
 * @author Adam Krejci
 */
class UniqueSequenceAlphabeticComparator implements Comparator<UniqueSequence>{

    @Override
    public int compare(UniqueSequence o1, UniqueSequence o2) {
        return o1.getSequenceString().compareTo(o2.getSequenceString());
    }   
}