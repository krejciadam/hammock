/*
 *Class represents result of extending clusters. Carries extended clusters
and sequences rejected because of violating minimum match states count
 */

package cz.krejciadam.hammock;

import java.util.List;

/**
 *
 * @author Adam Krejci
 */
public class ExtendedClusters {
    private final List<Cluster> clusters;
    private final List<UniqueSequence> rejectedSequences;

    public ExtendedClusters(List<Cluster> clusters, List<UniqueSequence> rejectedSequences) {
        this.clusters = clusters;
        this.rejectedSequences = rejectedSequences;
    }

    public List<Cluster> getClusters() {
        return clusters;
    }

    public List<UniqueSequence> getRejectedSequences() {
        return rejectedSequences;
    }
}
