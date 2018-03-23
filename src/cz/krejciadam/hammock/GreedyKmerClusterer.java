/*
 * Clusters sequences on the basis of shared kmers
 */
package cz.krejciadam.hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Adam Krejci
 */
public class GreedyKmerClusterer implements SequenceClusterer {

    KmerScorer scorer;
    int sharedKmers;
    int coreClustersCount;

    public GreedyKmerClusterer(KmerScorer scorer, int sharedKmers, int coreClustersCount) {
        this.scorer = scorer;
        this.sharedKmers = sharedKmers;
        this.coreClustersCount = coreClustersCount;
    }
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences) throws DataException{
        List<UniqueSequence> coreSequences = new ArrayList<>();
        List<UniqueSequence> remainingSequences = new ArrayList<>();
        for(UniqueSequence seq : sequences){
            int sum = 0;
            for (UniqueSequence core : coreSequences){
                Set<List<Integer>> coreKmers = scorer.getAllKmers(core);
                Set<List<Integer>> seqKmers = scorer.getAllKmers(seq);
                coreKmers.retainAll(seqKmers);
                sum += coreKmers.size();
            }
            if (sum == 0){
                coreSequences.add(seq);
            } else{
                remainingSequences.add(seq);
            }
        }
        List<Cluster> cores = new ArrayList<>();
        int id = 0;
        for(UniqueSequence seq : coreSequences){
            List<UniqueSequence> list = new ArrayList<>();
            list.add(seq);
            cores.add(new Cluster(list, id));
            id++;
        }
        List<Cluster> result = assignToClosestClusterIterative(cores, remainingSequences);
        for (Cluster cl : result){
            System.out.println(cl.getUniqueSize());
        }
        return(result);
    }


    
    private List<Cluster> assignToClosestClusterIterative(Collection<Cluster> clusters, Collection<UniqueSequence> sequences) throws DataException{
        List<Map<List<Integer>, Double>> kmerFrequencies = new ArrayList<>();
        List<Cluster> clusterList = new ArrayList<>(clusters);
        for (Cluster cl : clusterList) {
            Map<List<Integer>, Double> frequencyMap = new HashMap<>();
            for (UniqueSequence seq : cl.getSequences()) {
                Set<List<Integer>> kmers = scorer.getAllKmers(seq);
                for (List<Integer> kmer : kmers) {
                    Double count = frequencyMap.get(kmer);
                    if (count == null) {
                        count = 0.0;
                    }
                    count += 1.0;
                    frequencyMap.put(kmer, count);
                }
            }
            kmerFrequencies.add(frequencyMap);
        }
        
        for (UniqueSequence seq : sequences){
            Double bestScore = -1.0;
            int bestIndex = 0; //if no kmers present anywhere, goes to the first cluster
            for (int i = 0; i < clusterList.size(); i++){
                Double sum = 0.0;
                for (List<Integer> kmer : scorer.getAllKmers(seq)){
                    Double score = kmerFrequencies.get(i).get(kmer);
                    if (score != null){
                        sum += score;
                    }
                }
                sum = sum / (clusterList.get(i).getUniqueSize() + 0.0);
                if (sum > bestScore){
                    bestScore = sum;
                    bestIndex = i;
                }
            }
            clusterList.get(bestIndex).insert(seq);
            for (List<Integer> kmer : scorer.getAllKmers(seq)){
                    Double count = kmerFrequencies.get(bestIndex).get(kmer);
                    if (count == null) {
                        count = 0.0;
                    }
                    count += 1.0;
                    kmerFrequencies.get(bestIndex).put(kmer, count);
            }
        }
        return(clusterList);
    }
}
