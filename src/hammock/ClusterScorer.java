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
public interface ClusterScorer {
    public int clusterScore(Cluster cl1, Cluster cl2) throws DataException;
    
}
