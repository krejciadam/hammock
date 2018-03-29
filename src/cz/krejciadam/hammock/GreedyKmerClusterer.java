/*
 * Clusters sequences on the basis of shared kmers
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

/**
 *
 * @author Adam Krejci
 */
public class GreedyKmerClusterer implements SequenceClusterer {

    KmerScorer scorer;
    int sharedKmers;
    int coreClustersCount;
    ExecutorCompletionService<SearchResult> resultPool;

    public GreedyKmerClusterer(KmerScorer scorer, ExecutorService threadPool) {
        resultPool = new ExecutorCompletionService<>(threadPool);
        this.scorer = scorer;
    }
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences) throws DataException, InterruptedException, ExecutionException{
        List<UniqueSequence> coreSequences = new ArrayList<>();
        List<UniqueSequence> remainingSequences = new ArrayList<>();
        Set<List<Integer>> coreKmers = new HashSet<>();
        for(UniqueSequence seq : sequences){
            Set<List<Integer>> intersection = new HashSet<>(coreKmers);
            intersection.retainAll(scorer.getAllKmers(seq));
            if (intersection.isEmpty()){
                coreSequences.add(seq);
                coreKmers.addAll(scorer.getAllKmers(seq));
            } else{
                remainingSequences.add(seq);
            }
        }
        System.out.println(coreSequences.size());
        List<Cluster> cores = new ArrayList<>();
        int id = 0;
        for(UniqueSequence seq : coreSequences){
            List<UniqueSequence> list = new ArrayList<>();
            list.add(seq);
            cores.add(new Cluster(list, id));
            id++;
        }
        List<Cluster> result = assignToClosestClusterIterative(cores, remainingSequences);
        return(result);
    }


    
    private List<Cluster> assignToClosestClusterIterative(Collection<Cluster> clusters, Collection<UniqueSequence> sequences) throws DataException, InterruptedException, ExecutionException{
        List<Map<List<Integer>, Integer>> kmerFrequencies = new ArrayList<>();
        List<Cluster> clusterList = new ArrayList<>(clusters);
        for (Cluster cl : clusterList) {
            Map<List<Integer>, Integer> frequencyMap = new HashMap<>();
            for (UniqueSequence seq : cl.getSequences()) {
                Set<List<Integer>> kmers = scorer.getAllKmers(seq);
                for (List<Integer> kmer : kmers) {
                    Integer count = frequencyMap.get(kmer);
                    if (count == null) {
                        count = 0;
                    }
                    count += 1;
                    frequencyMap.put(kmer, count);
                }
            }
            kmerFrequencies.add(frequencyMap);
        }              
        for (UniqueSequence seq : sequences){
            for (int i = 0; i < clusterList.size(); i++){
                resultPool.submit(new bestClusterRunner(seq, i, kmerFrequencies.get(i), clusterList.get(i).getUniqueSize(), scorer));
            }
            Double bestScore = 0.0;
            int bestIndex = 0;
            for (int i = 0; i < clusterList.size(); i++){
                SearchResult res = resultPool.take().get();
                if (res.getScore() > bestScore){
                    bestScore = res.getScore();
                    bestIndex = res.getIndex();
                }
            }
            clusterList.get(bestIndex).insert(seq);
                        for (List<Integer> kmer : scorer.getAllKmers(seq)){
                    Integer count = kmerFrequencies.get(bestIndex).get(kmer);
                    if (count == null) {
                        count = 0;
                    }
                    count += 1;
                    kmerFrequencies.get(bestIndex).put(kmer, count);
            }
            
        }
        return(clusterList);
    }
}

class bestClusterRunner implements Callable<SearchResult>{
    
    private final UniqueSequence seq;
    private final int index;
    private final Map<List<Integer>, Integer> kmerFrequencies;
    private final int size;
    private final KmerScorer scorer;

    public bestClusterRunner(UniqueSequence seq, int index, Map<List<Integer>, Integer> kmerFrequencies, int size, KmerScorer scorer) {
        this.seq = seq;
        this.index = index;
        this.kmerFrequencies = kmerFrequencies;
        this.size = size;
        this.scorer = scorer;
    }
    
    @Override
    public SearchResult call() throws Exception {
        Double sum = 0.0;
        for (List<Integer> kmer : scorer.getAllKmers(seq)){
            Integer score = kmerFrequencies.get(kmer);
                if (score != null){
                    sum += score;
                }
            }
        sum = (sum + 0.0) / (size + 0.0);
        return(new SearchResult(index, sum));
    }
    
}

class SearchResult{
    private final int index;
    private final double score;

    public SearchResult(int index, double score) {
        this.index = index;
        this.score = score;
    }

    public int getIndex() {
        return index;
    }

    public double getScore() {
        return score;
    }
}
