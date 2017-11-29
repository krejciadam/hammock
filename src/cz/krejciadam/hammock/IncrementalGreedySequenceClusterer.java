/*
 * Class performs greedy clustering. Also aligns sequences to the first based on shift
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class IncrementalGreedySequenceClusterer implements SequenceClusterer {

    private final int threshold;

    /**
     * Constructor.
     *
     * @param threshold Minimum score for a sequence to be added into a cluster
     * during greedy clustering.
     */
    public IncrementalGreedySequenceClusterer(int threshold) {
        this.threshold = threshold;
    }

    /**
     * Performs greedy sequence clustering. 
     *
     * This method works in O(n^2) and is deterministic for defined input order of the sequences in sortedList
     *
     * @param sortedList Sequences to be sorted 
     * @param scorer scorer to be used
     * @return 
     * @throws InterruptedException
     * @throws ExecutionException
     */
    
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sortedList, AligningSequenceScorer scorer) throws InterruptedException, ExecutionException {
        int nThreads;
        CompletionService<AligningScorerResult> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);

        Map<UniqueSequence, List<UniqueSequence>> representativeMap = new HashMap<>(); //add first representative
        representativeMap.put(sortedList.get(0), new ArrayList<UniqueSequence>());
        List<UniqueSequence> representativeList = new ArrayList<>(); //defined outside map, for deterministic order (Map.keySet has artificial order and same-score may happen)
        representativeList.add(sortedList.get(0));
        
        for (int i = 1; i < sortedList.size(); i++) {
            nThreads = Hammock.nThreads;
            if (nThreads >= representativeList.size() + 2) {
                nThreads = 1;
            }
            UniqueSequence comparedSequence = sortedList.get(i);
            int forOneThread = representativeList.size() / nThreads; //*div* number division
            int j = 1;
            for (; j < nThreads; j++) {
                List<UniqueSequence> subList = representativeList.subList((j - 1) * forOneThread, j * forOneThread);
                resultPool.submit(new AligningBestHitRunner(subList, comparedSequence, scorer));
            }
            List<UniqueSequence> subList = representativeList.subList((j - 1) * forOneThread, representativeList.size()); //last piece of list
            resultPool.submit(new AligningBestHitRunner(subList, comparedSequence, scorer));
            AligningScorerResult bestResult = null;
            for (int t = 0; t < nThreads; t++) {
                AligningScorerResult partResult = resultPool.take().get();
                if ((bestResult == null) || (partResult != null && partResult.compareTo(bestResult) > 0)) {
                    bestResult = partResult;
                }
            }
            if (bestResult.getScore() >= threshold) {
                List<UniqueSequence> bestCluster = representativeMap.get(bestResult.getSequence());
                bestCluster.add(comparedSequence);
                representativeMap.put(bestResult.getSequence(), bestCluster);
            } else {
                representativeMap.put(comparedSequence, new ArrayList<UniqueSequence>());
                representativeList.add(comparedSequence);
            }
        }

        Map<Integer, UniqueSequence> pivotMap = new HashMap<>();
        Map<Integer, List<AligningScorerResult>> foundSeqsMap = new HashMap<>();
        List<Cluster> result = new ArrayList<>();
        int id = 1;

        
        for (UniqueSequence repreSeq : representativeList) {
            List<UniqueSequence> clusterSequences = new ArrayList<>();
            clusterSequences.add(repreSeq);
            clusterSequences.addAll(representativeMap.get(repreSeq));
            Cluster cl = new Cluster(clusterSequences, id);
            result.add(cl);
            pivotMap.put(id, repreSeq);
            List<AligningScorerResult> foundSeqsResult = new ArrayList<>();
            for (UniqueSequence seq : representativeMap.get(repreSeq)) {
                try{
                foundSeqsResult.add(scorer.scoreWithShift(repreSeq, seq));
                } catch(DataException e){
                    System.err.println("Error. This should not have happend.");
                }
            }
            foundSeqsMap.put(id, foundSeqsResult);
            id++;
        }
        Hammock.pivotMap = pivotMap;
        Hammock.foundSeqsMap = foundSeqsMap;
        return result;
    }
}

/**
 * Helper class for parallel processing. When searching through a list of
 * sequences with a pivot sequences, we divide the list in two parts - sequences
 * similar enough to pivot (selected) and the other sequences.
 *
 * @author Adam Krejci
 */
class AlignedListSearchResult implements Comparable<AlignedListSearchResult> {

    private final List<AligningScorerResult> selectedSequences;
    private final List<UniqueSequence> otherSequences;
    private final int order;

    public AlignedListSearchResult(List<AligningScorerResult> selectedSequences, List<UniqueSequence> otherSequences, int order) {
        this.selectedSequences = selectedSequences;
        this.otherSequences = otherSequences;
        this.order = order; //we remember the order, so we do not need to sort.
    }

    public List<AligningScorerResult> getSelectedSequences() {
        return selectedSequences;
    }

    public List<UniqueSequence> getOtherSequences() {
        return otherSequences;
    }
 
    public int getOrder() {
        return order;
    }

    @Override
    public int compareTo(AlignedListSearchResult o) {
        return this.getOrder() - o.getOrder();
    }
}

/**
 * Helper class for parallel processing. Searches through a given list of
 * sequences using a pivot sequences. All sequences similar to pivot above given
 * threshild are reported as selected, all the other sequences are reported as
 * other sequences.
 *
 * @author Adam Krejci
 */
class AlignedListSearchRunner implements Callable<AlignedListSearchResult> {

    private final List<UniqueSequence> sequenceList;
    private final UniqueSequence pivot;
    private final int scoreThreshold;
    private final AligningSequenceScorer scorer;
    private final int order; //pamatujeme si poradi, abychom zabranili nutnosti trideni

    public AlignedListSearchRunner(List<UniqueSequence> phageList, UniqueSequence pivot, int scoreThreshold, AligningSequenceScorer scorer, int order) {
        this.sequenceList = phageList;
        this.pivot = pivot;
        this.scoreThreshold = scoreThreshold;
        this.scorer = scorer;
        this.order = order;
    }

    @Override
    public AlignedListSearchResult call() throws DataException {
        AligningScorerResult result;
        List<AligningScorerResult> selectedSequences = new ArrayList<>();
        List<UniqueSequence> otherSequences = new ArrayList<>();
        for (UniqueSequence comparedSequence : sequenceList) {
            result = scorer.scoreWithShift(pivot, comparedSequence);
            if (result.getScore() >= scoreThreshold) {
                selectedSequences.add(result);
            } else {
                otherSequences.add(comparedSequence);
            }
        }
        return new AlignedListSearchResult(selectedSequences, otherSequences, order);
    }
}


/**
 * Helper class for parallel processing. Compares comparedSequence to all sequences
 * in databaseList and returns the best hit.
 * @author Adam Krejci
 */
class AligningBestHitRunner implements Callable<AligningScorerResult> {

    private final List<UniqueSequence> databaseList;
    private final UniqueSequence comparedSequence;
    private final AligningSequenceScorer scorer;

    public AligningBestHitRunner(List<UniqueSequence> databaseList, UniqueSequence comparedSequence, AligningSequenceScorer scorer) {
        this.databaseList = databaseList;
        this.comparedSequence = comparedSequence;
        this.scorer = scorer;
    }

    @Override
    public AligningScorerResult call() throws Exception {
        AligningScorerResult bestResult = null;
        for (UniqueSequence seq : databaseList) {
            AligningScorerResult roundResult = scorer.scoreWithShift(comparedSequence, seq);
            if (bestResult == null || (roundResult.compareTo(bestResult) > 0)) {
                bestResult = roundResult;
            }
        }
        return bestResult;
    }
}
