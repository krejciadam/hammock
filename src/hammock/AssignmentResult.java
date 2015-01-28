/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package hammock;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Class represents results of one clustering iteration - clusters and remaining
 * (orphan) database sequences.
 *
 * @author Adam Krejci
 */
class AssignmentResult {

    private final List<Cluster> clusters = new ArrayList<>();
    private final List<UniqueSequence> databaseSequences = new ArrayList<>();

    public AssignmentResult(Collection<Cluster> clusters, Collection<UniqueSequence> databaseSequences) {
        this.clusters.addAll(clusters);
        this.databaseSequences.addAll(databaseSequences);
    }

    public List<Cluster> getClusters() {
        return clusters;
    }

    public List<UniqueSequence> getDatabaseSequences() {
        return databaseSequences;
    }
}