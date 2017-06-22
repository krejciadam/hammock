/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;

/**
 *
 * @author Adam Krejci
 */
public class ClinkageSequenceClusterer implements SequenceClusterer {


    private final int sizeLimit; //do pameti budeme ukladat skore jen pro takto velke
                                //a vetsi clustery (jen pro CachedCLScorer)
    private final int threshold;

    public ClinkageSequenceClusterer(int threshold) {
        this.threshold = threshold;
        this.sizeLimit = 1;
    }
    
    public ClinkageSequenceClusterer(int threshold, int sizeLimit) {
        this.threshold = threshold;
        this.sizeLimit = sizeLimit;
    }
     
    
    @Override
    public List<Cluster> cluster(List<UniqueSequence> sequences, AligningSequenceScorer sequenceScorer) throws InterruptedException, ExecutionException {
        CachedClusterScorer clusterScorer = new CachedClusterScorer(new ClinkageClusterScorer(sequenceScorer, threshold), sizeLimit);
        ClusterStack<Cluster> stack = new ClusterStack<>();
        Set<Cluster> readyClusters = new HashSet<>();  //clustery, jejichz nejblizsi soused je za thresholdem budeme vyndavat do tohoto listu
        int sumCommodity = 0; //Celkova cena vsech zbyvajiich clusteru    
        Set<Cluster> activeClusters = new HashSet<>();
        int currentId = 1;
        for (UniqueSequence seq : sequences){
            Set<UniqueSequence> seqs = new HashSet<>();
            seqs.add(seq);
            activeClusters.add(new Cluster(seqs, currentId));
            currentId ++;
        }
        
        /* Inicialni zjisteni celkove ceny vsech clusteru */
        for (Cluster cluster : activeClusters) {
            sumCommodity += cluster.getUniqueSize();
        }
        
        CompletionService<NearestCluster> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        
        while (activeClusters.size() > 1) {
            System.out.print("\rSequences remaining to process: "
                    + activeClusters.size() + ", clusters ready: "
                    + readyClusters.size() + " ");

            Cluster randomCluster = activeClusters.iterator().next();
            stack.push(randomCluster); // pridame libovolny prvek
            while (!stack.isEmpty()) {
                Cluster top = stack.peek();
                
                /*Hledani nejblizsiho clusteu*/
                
                NearestCluster foundNearestCluster = findNearestClusterParallel(activeClusters, top, clusterScorer, sumCommodity, resultPool);
                int maxScore = Integer.MIN_VALUE;
                Cluster nearestCluster = null;
                if (foundNearestCluster != null){
                    nearestCluster = foundNearestCluster.getCluster();
                    maxScore = foundNearestCluster.getScore();
                }
                
                 /*pokud je skore moc male, odstranime top a jedeme od znova (toto muze nastat jen kdyz je v zasobniku pouze jeden cluster)*/
                if (maxScore < threshold) {
                    stack.pop();
                    readyClusters.add(top);
                    activeClusters.remove(top);
                    sumCommodity -= top.getUniqueSize(); //snizeni celkove ceny aktivnich clusteru
                    continue;
                }
                
                /*Pokud je cluster nejblizsi topu shodny s 2. clusterem ve stacku
                 * Vyjmeme oba ze stacku a z aktivnich clusteru, sloucime je a vlozime
                 slouceny cluster do aktivnich clusteru. activeInputs se nemeni*/
                if ((stack.size() > 1) && stack.peekSecond().equals(nearestCluster)) { //prima rovnost: muze byt equals?
                    currentId++;
                    stack.pop();
                    stack.pop();
                    activeClusters.remove(top);
                    activeClusters.remove(nearestCluster);
                    clusterScorer.join(top, nearestCluster, currentId);
                    sumCommodity -= (top.getUniqueSize() + nearestCluster.getUniqueSize()); //snizime celkovou cenu o odebrane clustery
                    
                    List<UniqueSequence> mergedSeqs = top.getSequences();
                    mergedSeqs.addAll(nearestCluster.getSequences());
                    
                    Cluster newTop = new Cluster(mergedSeqs, currentId);

                    activeClusters.add(newTop);
                    sumCommodity += newTop.getUniqueSize(); //zvysime celkovou cenu o nove pridany cluster //??? (-= bylo spatne?)
                } else { //ve stacku je jen 1 prvek nebo druhy prvek neni roven prvnimu
                    stack.push(nearestCluster);
                }
            }
        }
        /*presuneme posledni cluster do hotovych*/
        readyClusters.add(activeClusters.iterator().next());

        /*vracime list*/
        List<Cluster> toReturn = new ArrayList<>();
        toReturn.addAll(readyClusters);
        return new ArrayList<>(toReturn);
    }
    
    public static NearestCluster findNearestClusterParallel(Collection<Cluster> inputClusters, Cluster comparedCluster, ClusterScorer scorer, int sumCommodity, CompletionService<NearestCluster> resultPool) throws InterruptedException, ExecutionException{
        if(inputClusters.isEmpty()){
            return(new NearestCluster(null, Integer.MIN_VALUE));
        }
        int currentMaxNumberOfThreads = ClinkageSequenceClusterer.setNumberOfParts(inputClusters);
        List<Set<Cluster>> activeClustersParts = ClinkageSequenceClusterer.activeClustersParts(inputClusters, sumCommodity, currentMaxNumberOfThreads);
        /*Vlastni hledani nejblizsiho clusteru*/
        for (int i = 0; i < activeClustersParts.size(); i++) {
            resultPool.submit(
                new NearestClusterRunner(
                activeClustersParts.get(i),
                comparedCluster, scorer));
            }
                
        int maxScore = (Integer.MIN_VALUE + 42); //kdyz je NearestCluster z dalsi casti NULL, ma MIN_VALUE. Takto zajistime, ze takove pripady se rovnou preskoci
        NearestCluster nearestCluster = null;
                
            /*Ziskani vysledku z jiz dopocitanych vlaken*/
        for (Set<Cluster> activeClustersPart : activeClustersParts) {
            NearestCluster currentCluster = resultPool.take().get();     
            
            
            if (currentCluster.getScore() < maxScore) {
                continue;
            }

            if (currentCluster.getScore() > maxScore) {
                nearestCluster = currentCluster;
                maxScore = nearestCluster.getScore();
            } else { //currentScore == maxScore   
                if ((currentCluster.getCluster().size() > nearestCluster.getCluster().size())) { //kdyz je current cluster null (vyjimecna situace, kdyz byl v 1 vlakne jen top), muzeme preskocit
                    nearestCluster = currentCluster;
                } else{
                    if ((currentCluster.getCluster().size() == nearestCluster.getCluster().size()) && (currentCluster.getCluster().getId() < nearestCluster.getCluster().getId())){ //shodna velikost
                        nearestCluster = currentCluster;
                    }
                }
            }
        }
        return(nearestCluster);
    }
    
    
        /**
     * Urcuje, na kolik casti budeme delit praci v zavislosti na aktualnim poctu
     * zpracovavanych clusteru.
     *
     * @param activeClusters Mnozina zpracovavanych clusteru
     * @return
     */
    private static int setNumberOfParts(Collection<Cluster> activeClusters) {
        int currentMaxNumberOfParts = Hammock.nThreads * 4;
        if (activeClusters.size() < Hammock.nThreads * 4 + 1) {          //Kdyz uz je malo aktivnich clusteru, musime delat i mene casti
            currentMaxNumberOfParts = Math.max(activeClusters.size() - 1, 1);
        }
        return currentMaxNumberOfParts;
    }
    
        /**
     * rozdeleni mnoziny aktivnich clusteru na casti (jedna cast pro kazde
     * vlakno) Snazime se o spravedlive rozdeleni, tedy aby kazde vlakno melo
     * nejlepe stejnou porci "komodity"
     *
     * @param activeClusters Clustery, ktere budeme rozdelovat na stejne tezke
     * casti
     * @param sumCommodity Celkova suma "komodity narocnosti" v activeClusters
     * @param numberOfParts pocet casti, na ktere budeme praci delit
     * @return Seznam mnozin clusteru. Kazda mnozina je porca prace pro 1 vlakno
     */
    private static List<Set<Cluster>> activeClustersParts(Collection<Cluster> activeClusters,
            long sumCommodity,
            int numberOfParts) {

        List<Set<Cluster>> result = new ArrayList<>();
        Set<Cluster> currentSet = new HashSet<>();
        long forOneThread = (sumCommodity / numberOfParts) + 1; //celociselne deleni
        long portionCounter = forOneThread;
        for (Cluster clust : activeClusters) {
            currentSet.add(clust);
            portionCounter -= clust.getUniqueSize();
            if (portionCounter <= 0) {
                result.add(currentSet);
                currentSet = new HashSet<>();
                portionCounter = forOneThread;
            }
        }
        if (currentSet.size() > 0) {
            result.add(currentSet); //pridame i posledni set clusteru
        }
        return result;
    }
    
    
}

class ClusterStack<E> extends Stack<E> {
    /**
     * Vrati element, ktery je v zasobniku druhy (prvni pod vrcholem)
     * @return druhy element zasobniku
     */
    public E peekSecond() {
        E element = this.get(this.size() - 2);
        return element;
    }
}


/**
 * Callable pro hledani nejpodobnejsiho clusteru pri vypoctu shlukovaciho
 * algoritmu.
 */
class NearestClusterRunner implements Callable<NearestCluster> {

    private final Collection<Cluster> databaseClusters;
    private final Cluster comparedCluster;
    private final ClusterScorer scorer;

    public NearestClusterRunner(Collection<Cluster> activeClustersPart,
            Cluster top,
            ClusterScorer scorer) {
        this.databaseClusters = activeClustersPart;
        this.comparedCluster = top;
        this.scorer = scorer;
    }

    @Override
    public NearestCluster call() throws DataException{
        int maxScore = Integer.MIN_VALUE;
        int score;
        Cluster nearestCluster = null;
        for (Cluster i : databaseClusters) {
            score = scorer.clusterScore(i, comparedCluster);
            if (score < maxScore) {
                continue;
            }
            if (score > maxScore) {
                if ((i != comparedCluster)) {
                    maxScore = score;
                    nearestCluster = i;
                }
            } /*else vetev: (score == maxScore) Pokud existuje vice nejblizsich clusteru,
             * musime se DETERMINISTICKY rozhodnout, ktery prohlasime za
             * nejblizsi. Bereme pak jednoznacne kriterium velikosti. Pokud i ta je stejna,
             * usporadani podle klice: atribut id, viz. algoritmus nearest neighbour chain.
             */ else {
                if (i != comparedCluster){
                    if (i.size() > nearestCluster.size()){
                        nearestCluster = i;                        
                    } else{
                        if (i.size() < nearestCluster.size()){
                           //nic - zadna zmena 
                        }
                        else { //i velikost je stejna, pak se rozhodneme podle id
                            if ((i.getId() < nearestCluster.getId())) {
                            nearestCluster = i;
                        }    
                        }
                    }
                }
            }
        }
        return new NearestCluster(nearestCluster, maxScore);
    }
}