/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hammock;

/**
 *
 * @author akrejci
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
