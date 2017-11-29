/*
 * A class to represent the results from a cluster-cluster comparison search.
 */
package cz.krejciadam.hammock;

/**
 *
 * @author Adam Krejci
 */
public class NearestCluster {

    private Cluster cluster;
    private int score;

    public NearestCluster(Cluster cluster, int score) {
        this.cluster = cluster;
        this.score = score;
    }

    public Cluster getCluster() {
        return cluster;
    }

    public void setCluster(Cluster cluster) {
        this.cluster = cluster;
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }
}
