/*
 * Class represents signle result of HHalign run - score of alignment
   of two hmms representing two clusters.
 */

package hammock;

import java.util.Objects;

/**
 *
 * @author Adam Krejci
 */
public class HHalignHit implements Comparable<HHalignHit> {
    
    private final UnorderedPair<Cluster> pair;
    private final Double score;
    private final Cluster searchedCluster;
    private final Cluster foundCluster;
    private final String alignmentLine1;
    private final String alignmentLine2;

    public HHalignHit(Cluster searchedCluster, Cluster foundCluster, Double score, String alignmentLine1, String alignmentLine2) {
        this.pair = new UnorderedPair<>(searchedCluster, foundCluster);
        this.score = score;
        this.searchedCluster = searchedCluster;
        this.foundCluster = foundCluster;
        this.alignmentLine1 = alignmentLine1;
        this.alignmentLine2 = alignmentLine2;
    }

    public Cluster getBiggerCluster() {
        return pair.getBigger();
    }

    public Cluster getSmallerCluster() {
        return pair.getSmaller();
    }

    public UnorderedPair<Cluster> getPair() {
        return pair;
    }
    

    public Double getScore() {
        return score;
    }

    public Cluster getSearchedCluster() {
        return searchedCluster;
    }

    public Cluster getFoundCluster() {
        return foundCluster;
    }

    public String getAlignmentLine1() {
        return alignmentLine1;
    }

    public String getAlignmentLine2() {
        return alignmentLine2;
    }
    
    
    
    

    @Override
    public int compareTo(HHalignHit o) {
        if (this.pair.equals(o.getPair())){
            return 0;
        }
        if (this.score > o.getScore()){
            return 1;
        } if (this.score < o.getScore()){
            return -1;
        }
        int sizeDifference = (this.getBiggerCluster().size() + this.getSmallerCluster().size()) - (o.getBiggerCluster().size() + o.getSmallerCluster().size());
        if (sizeDifference != 0){
            return sizeDifference;
        } else{
            return this.getBiggerCluster().getSequences().get(0).getSequenceString().compareTo(o.getBiggerCluster().getSequences().get(0).getSequenceString());
        }   
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 41 * hash + Objects.hashCode(this.pair);
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
        final HHalignHit other = (HHalignHit) obj;
        if (!Objects.equals(this.pair, other.pair)) {
            return false;
        }
        return true;
    }

    
}
