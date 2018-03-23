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
    
    public List<Cluster> cluster2(List<UniqueSequence> sequences) throws DataException{
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
        for (Cluster cl : result){
            System.out.println(cl.getUniqueSize());
        }
        return(result);
    }


    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences) throws InterruptedException, ExecutionException, DataException {
        List<Cluster> res = new ArrayList<>();
        int id = 0;
        Set<UniqueSequence> remainingSequences = new HashSet<>();
        remainingSequences.addAll(sequences);

        for (int j = 0; j < coreClustersCount; j++) {
            Map<List<Integer>, Set<UniqueSequence>> kmerSeqMap = makeKmerSeqMap(remainingSequences);
            Set<List<Integer>> usedKmers = new HashSet<>();
            Set<UniqueSequence> clusterSequences = new HashSet<>();
            for (int i = 0; i < sharedKmers; i++) {
                List<Integer> maxKmer = null;
                int maxCount = 0;
                for (Map.Entry<List<Integer>, Set<UniqueSequence>> entry : kmerSeqMap.entrySet()) {
                    if ((entry.getValue().size() > maxCount) && !(usedKmers.contains(entry.getKey()))) { //has not been counted already
                        maxKmer = entry.getKey();
                        maxCount = entry.getValue().size();
                    }
                }
                usedKmers.add(maxKmer);
                clusterSequences = kmerSeqMap.get(maxKmer);
                if (i < sharedKmers - 1) {
                    kmerSeqMap = makeKmerSeqMap(kmerSeqMap.get(maxKmer));
                }
            }
            res.add(new Cluster(clusterSequences, id));
            id++;
            remainingSequences.removeAll(clusterSequences);
        }
        for (UniqueSequence seq : remainingSequences) {
            List<UniqueSequence> seqList = new ArrayList<>();
            seqList.add(seq);
            res.add(new Cluster(seqList, id));
            id++;
        }
        return (res);
    }
    
    private List<Cluster> assignToClosestClusterIterative(Collection<Cluster> clusters, Collection<UniqueSequence> sequences) throws DataException{
        List<Map<List<Integer>, Double>> kmerFrequencies = new ArrayList<>();
        List<Cluster> clusterList = new ArrayList<>(clusters);
        for (Cluster cl : clusterList) {
            Map<List<Integer>, Double> frequenceMap = new HashMap<>();
            for (UniqueSequence seq : cl.getSequences()) {
                Set<List<Integer>> kmers = scorer.getAllKmers(seq);
                for (List<Integer> kmer : kmers) {
                    Double count = frequenceMap.get(kmer);
                    if (count == null) {
                        count = 0.0;
                    }
                    count += 1.0;
                    frequenceMap.put(kmer, count);
                }
            }
            kmerFrequencies.add(frequenceMap);
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

    public List<Cluster> assignToClosestCluster(Collection<Cluster> clusters, Collection<UniqueSequence> sequences) throws DataException {
        List<Map<List<Integer>, Double>> kmerFrequencies = new ArrayList<>();
        List<Cluster> clusterList = new ArrayList<>(clusters);
        for (Cluster cl : clusterList) {
            Map<List<Integer>, Double> frequenceMap = new HashMap<>();
            for (UniqueSequence seq : cl.getSequences()) {
                Set<List<Integer>> kmers = scorer.getAllKmers(seq);
                for (List<Integer> kmer : kmers) {
                    Double count = frequenceMap.get(kmer);
                    if (count == null) {
                        count = 0.0;
                    }
                    count += 1.0;
                    frequenceMap.put(kmer, count);
                }
            }
            for (Map.Entry<List<Integer>, Double> entry : frequenceMap.entrySet()) {
                frequenceMap.put(entry.getKey(), entry.getValue() / cl.getUniqueSize());
            }
            kmerFrequencies.add(frequenceMap);
        }
        List<List<UniqueSequence>> toAdd = new ArrayList<>();
        for (int i = 0; i < clusterList.size(); i++){
            toAdd.add(new ArrayList<UniqueSequence>());
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
                if (sum > bestScore){
                    bestScore = sum;
                    bestIndex = i;
                }
            }
            toAdd.get(bestIndex).add(seq);
        }
        for (int i = 0; i < clusterList.size(); i++){
            clusterList.get(i).insertAll(toAdd.get(i));
        }
        return(clusterList);
    }

    private Map<List<Integer>, Set<UniqueSequence>> makeKmerSeqMap(Collection<UniqueSequence> sequences) {
        Map<List<Integer>, Set<UniqueSequence>> kmerSeqMap = new HashMap<>();
        for (UniqueSequence seq : sequences) {
            Set<List<Integer>> kmers = scorer.getAllKmers(seq);
            for (List<Integer> kmer : kmers) {
                Set<UniqueSequence> seqSet = kmerSeqMap.get(kmer);
                if (seqSet == null) {
                    seqSet = new HashSet<>();
                }
                seqSet.add(seq);
                kmerSeqMap.put(kmer, seqSet);
            }
        }
        return (kmerSeqMap);
    }

}
