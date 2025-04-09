/*
 * Class performs Hmm-Hmm comparison tasks by running HHSuite as external process 
 */
package cz.krejciadam.hammock;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

/**
 *
 * @author Adam Krejci
 */
public class HHsuiteRunner {

    /**
     * Builds hmms in hhsuite format for collection of clusters using parallel
     * execution of multiple instances of hhmake. Clusters with hasHH == ture
     * will not have their HMMs built.
     *
     * Sets hasHH to true to every cluster if successful. Clusters must be
     * modifiable. Writes out files with hmms to hh directory specified by
     * Settings. File name is in format: [cluster_id].hhm If such file exists
     * for a cluster having hasHH() set to false, file will be overwritten by a
     * new file without warning.
     *
     * @param clusters clusters for which hmms will be built.
     * @param threadPool ExecutorService to be used for parallel run
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static void buildHHs(Collection<Cluster> clusters, ExecutorService threadPool) throws InterruptedException, ExecutionException, Exception {
        List<Callable<Void>> callables = new ArrayList<>();
        for (Cluster cluster : clusters) {
            callables.add(new SingleThreadHHmakeRunner(cluster));
        }
        ParallelRunner.ExecuteInParallel(callables, threadPool);
    }

    /**
     * Builds hmm in hhsuite format for single cluster. If the cluster already
     * hasHH(), nothing happens. Sets hasHH to true to this cluster if
     * successful. Cluster must be modifiable. Writes out file with hmms to hh
     * directory specified by Settings. File name is in format: [cluster_id].hhm
     * If such file exists for a cluster having hasHH() set to false, file will
     * be overwritten by a new file with no warning.
     *
     * @param cluster cluster for which hmm will be built
     * @throws Exception
     */
    public static void buildHH(Cluster cluster) throws Exception {
        Callable callable = new SingleThreadHHmakeRunner(cluster);
        callable.call();
    }

    public static int getHHLength(Cluster cl) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(
                Settings.getInstance().getHhDirectory() + cl.getId() + ".hhm")))) {
            String line = null;
            for (int lineNum = 0; lineNum < 7; lineNum++) {
                line = reader.readLine();
            }
            String[] splitLine = line.split("\\s+"); //split on any number of whitespace characters
            return Integer.decode(splitLine[1]);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Compares alignedCluster to all members of searchedClusters. Returns all
     * identified hits (possibly less hits than size of searchedClusters)
     *
     * @param alignedCluster Cluster to be compared to all members of
     * searchedClusters
     * @param searchedClusters A collection of clusters to be compared to alignedCluster
     * @paragetHHLengtm searchedClusters All members of this Collection will be compared
     * to alignedCluster
     * @param threadPool ExecutorService to be used for parallel run
     * @return list of all identified hits. This list may be smaller than the
     * size of searchedClusters (hits with poor score are not always reported)
     * @throws Exception
     */
    public static List<HHalignHit> alignHmmList(Cluster alignedCluster, Collection<Cluster> searchedClusters, ExecutorService threadPool) throws Exception {
        List<List<Cluster>> superList = new ArrayList<>();
        int listCount = Math.min(Hammock.nThreads, searchedClusters.size());
        for (int i = 0; i < listCount; i++) {
            superList.add(new ArrayList<Cluster>());
        }
        int clustersAdded = 0;
        for (Cluster cl : searchedClusters) {
            superList.get(clustersAdded % listCount).add(cl);
            clustersAdded++;
        }
        CompletionService<List<HHalignHit>> resultPool = new ExecutorCompletionService<>(threadPool);
        for (int i = 0; i < superList.size(); i++) {
            resultPool.submit(new SingleThreadHHsearchRunner(alignedCluster, superList.get(i), "" + alignedCluster.getId() + "_" + i));
        }
        List<HHalignHit> result = new ArrayList<>();
        for (int i = 0; i < superList.size(); i++) {
            result.addAll(resultPool.take().get());
        }
        return result;
    }
    
    public static HHalignHit alignClusters(Cluster alignedCluster, Cluster searchedCluster) throws Exception{
        List<Cluster> tempList = new ArrayList<>();
        tempList.add(searchedCluster);
        Callable<List<HHalignHit>> callable = new SingleThreadHHsearchRunner(alignedCluster, tempList, "" + alignedCluster.getId() + "_" + 0);
        List<HHalignHit> resList = callable.call();
        return(resList.get(0));
    }

    /**
     * Performs all vs. all comparison of clusters using hhserach. Returns all
     * identified hits.
     *
     * @param clusters Clusters to be compared to each other
     * @param threadPool ExecutorService to be used for parallel run
     * @return All identifed hits
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static List<HHalignHit> alignAllVsAll(Collection<Cluster> clusters, ExecutorService threadPool) throws InterruptedException, ExecutionException {
        List<Cluster> clusterList = new ArrayList<>();
        clusterList.addAll(clusters);
        CompletionService<List<HHalignHit>> resultPool = new ExecutorCompletionService<>(threadPool);
        for (int i = 0; i < clusterList.size(); i++) {
            resultPool.submit(new SingleThreadHHsearchRunner(clusterList.get(i), clusterList.subList(i + 1, clusterList.size()), "" + clusterList.get(i).getId()));
        }
        List<HHalignHit> result = new ArrayList<>();
        for (int i = 0; i < clusterList.size(); i++) {
            result.addAll(resultPool.take().get());
        }
        return result;
    }
    
    /**
     * Aligns all clusters from clusters1 to all clusters from clusters2
     * @param clusters1 One collection of clusters to align
     * @param clusters2 Another collection of clusters to align
     * @param threadPool To be used for parallel prcessing
     * @return All the identified hits
     * @throws InterruptedException
     * @throws ExecutionException 
     */
    public static List<HHalignHit> alignAllVsAll(Collection<Cluster> clusters1, Collection<Cluster> clusters2, ExecutorService threadPool) throws InterruptedException, ExecutionException{
        List<Cluster> clusterList = new ArrayList<>();
        clusterList.addAll(clusters2);
        CompletionService<List<HHalignHit>> resultPool = new ExecutorCompletionService<>(threadPool);
        for (Cluster cl: clusters1){
            resultPool.submit(new SingleThreadHHsearchRunner(cl, clusterList, "" + cl.getId()));
        }
        List<HHalignHit> result = new ArrayList<>();
        for (Cluster cl : clusters1){
            result.addAll(resultPool.take().get());
        }
        return result;
    }
    
    /**
     * Merges two clusters contained in a HHalignHit object and returns resulting
     * cluster having a MSA constructed. Method uses alignment supported by
     * HHalign/HHsearch
     *
     * @param hit object representing two clusters to be merged
     * @param newId Resulting merged cluster will have this id
     * @return Cluster with MSA containing all sequences from both clusters on
     * input.
     * @throws DataException
     * @throws IOException
     */
    public static Cluster mergeClusters(HHalignHit hit, int newId) throws DataException, IOException {
        return mergeClusters(hit.getSearchedCluster(), hit.getFoundCluster(), hit.getAlignmentLine1(), hit.getAlignmentLine2(), newId);
    }

    /**
     * Merges two clusters and returns resulting cluster having MSA. Appropriate
     * lines from file containing result of HHsearch/HHalign run are required.
     *
     * @param cl1 First cluster to be merged
     * @param cl2 second cluster to be merged
     * @param resultLine1 first alignment line from HHalign/HHserach result
     * @param resultLine2 second alignment line form HHalign/HHserach result
     * @param newId Resulting cluster will have this id
     * @return Cluster with MSA containing all sequences from both clusters on
     * input.
     * @throws DataException
     * @throws IOException
     */
    private static Cluster mergeClusters(Cluster cl1, Cluster cl2, String resultLine1, String resultLine2, int newId) throws DataException, IOException {
        List<List<Integer>> gaps = getNewGapPositions(resultLine1, resultLine2, FileIOManager.getFirstA2mSeq(cl1), FileIOManager.getFirstA2mSeq(cl2));
        FileIOManager.mergeAlignedClusters(cl1, cl2, gaps.get(0), gaps.get(1), newId);
        Cluster newCluster = new Cluster(cl1.getSequences(), newId);
        newCluster.insertAll(cl2.getSequences());
        newCluster.setAsHasMSA();
        return newCluster;
    }

    /**
     * Parses HHalign/HHsearch result and finds out where to add gaps to MSA's
     * representing two clusters aligned.
     *
     * @param resultLine1 first alignment line from HHalign/HHserach result
     * @param resultLine2 second alignment line form HHalign/HHserach result
     * @param alignedSeq1 first line of MSA of first cluster
     * @param alignedSeq2 first line of MSA of second cluster
     * @return list containing two lists, first contains positions of new gaps
     * to be inserted into first cluster's MSA, second the same for second
     * cluster.
     */
    private static List<List<Integer>> getNewGapPositions(String resultLine1, String resultLine2, String alignedSeq1, String alignedSeq2) {
        String[] splitLine1 = resultLine1.split("\\s+"); //split on any number of whitespace characters
        String[] splitLine2 = resultLine2.split("\\s+");
        int startingLetter1 = Integer.parseInt(splitLine1[2]);
        int startingLetter2 = Integer.parseInt(splitLine2[2]);
        String alignmentLine1 = splitLine1[3];
        String alignmentLine2 = splitLine2[3];
        int endingLetter1 = Integer.parseInt(splitLine1[4]);
        int endingLetter2 = Integer.parseInt(splitLine2[4]);

        List<Integer> newGaps1 = new ArrayList<>();
        List<Integer> newGaps2 = new ArrayList<>();

        int lettersCounted1 = 0;
        int position1 = 0;
        while (lettersCounted1 < startingLetter1) {
            if (alignedSeq1.charAt(position1) != '.' && alignedSeq1.charAt(position1) != '-') {
                lettersCounted1++;
            }
            position1++;
            if (lettersCounted1 == endingLetter1){ // avoids null pointer exception when resultLine is all gaps
                break;
            }
        }

        int lettersCounted2 = 0;
        int position2 = 0;
        while (lettersCounted2 < startingLetter2) {
            if (alignedSeq2.charAt(position2) != '.' && alignedSeq2.charAt(position2) != '-') {
                lettersCounted2++;
            }
            position2++;
            if (lettersCounted2 == endingLetter2){ // avoids null pointer exception when resultLine is all gaps
                break;
            }
        }

        if (position1 != position2) { //need to add gaps to beginining
            for (int i = 0; i < Math.abs(position2 - position1); i++) {
                if (position1 < position2) {
                    newGaps1.add(i);
                } else {
                    newGaps2.add(i);
                }
            }
        }

        position1--;
        position2--;
        int offset1 = newGaps1.size();
        int offset2 = newGaps2.size();

        //add inside gaps if any
        for (int i = 0; i < alignmentLine1.length(); i++) {
            if (alignedSeq1.length() > (position1 + i)) {
                if ((alignmentLine1.charAt(i) == '-') && (alignedSeq1.charAt(position1 + i) != '.')) {
                    newGaps1.add(position1 + i + offset1);
                }
            } else {
                newGaps1.add(position1 + i + offset1);
            }
        }

        for (int i = 0; i < alignmentLine2.length(); i++) {
            if (alignedSeq2.length() > (position2 + i)) {
                if ((alignmentLine2.charAt(i) == '-') && (alignedSeq2.charAt(position2 + i) != '.')) {
                    newGaps2.add(position2 + i + offset2);
                }
            } else {
                newGaps2.add(position2 + i + offset2);
            }
        }

        // add gaps to end
        int seqLength1 = alignedSeq1.length() + newGaps1.size();
        int seqLength2 = alignedSeq2.length() + newGaps2.size();

        if (seqLength1 < seqLength2) {
            for (int i = seqLength1; i < seqLength2; i++) {
                newGaps1.add(i);
            }
        }

        if (seqLength2 < seqLength1) {
            for (int i = seqLength2; i < seqLength1; i++) {
                newGaps2.add(i);
            }
        }

        List<List<Integer>> result = new ArrayList<>();
        result.add(newGaps1);
        result.add(newGaps2);
        return result;
    }
}

/**
 * Runs hhmake on sigle thread as external process.
 *
 * @author Adam Krejci
 */
class SingleThreadHHmakeRunner implements Callable<Void> {

    private final Cluster cluster;

    /**
     * @param cluster Cluster for which hhsuite hmm will be built
     */
    public SingleThreadHHmakeRunner(Cluster cluster) {
        this.cluster = cluster;
    }

    /**
     * Builds a HMM in hhsuite format for single cluster. Runs hhmake as
     * external process. Resulting file will be saved in hh folder specified by
     * Settings. File name format: [cluster_id].hhm If cluster hasHH(), nothing
     * happens. If cluster.hasMSA() == false, MSA is constructed first. Sets
     * cluster's hasHH() to if successful.
     *
     * cluster must be modifiable.
     *
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    @Override
    public Void call() throws IOException, InterruptedException, DataException {
        List<String> parameters = new ArrayList<>();
        if (cluster.hasHH()) {
            return null;
        }
        if (!(cluster.hasMSA())) {
            ClustalRunner.multipleAlignment(cluster);
        }
        FileIOManager.alnToA2M(cluster);
        parameters.add("-i");
        parameters.add(Settings.getInstance().getMsaDirectory() + cluster.getId() + ".a2m");
        parameters.add("-o");
        parameters.add(Settings.getInstance().getHhDirectory() + cluster.getId() + ".hhm");
        parameters.add("-name");
        parameters.add(Integer.toString(cluster.getId()));
        if (Settings.getInstance().getHhmakeParameters() != null) {
            parameters.addAll(Settings.getInstance().getHhmakeParameters());
        }
        ExternalProcessRunner.runProcess(Settings.getInstance().getHhmakeCommand(), parameters, null, System.err, "While builidng HH for cluster " + cluster.getId() + " ", Hammock.hhSuiteEnv);
        cluster.setAsHasHH();
        return null;
    }
}

/**
 * Runs hhsearch in single thread as external process
 *
 * @author Adam Krejci
 */
class SingleThreadHHsearchRunner implements Callable<List<HHalignHit>> {

    private final Cluster alignedCluster;
    private final Collection<Cluster> searchedClusters;
    private final String idString;

    /**
     *
     * @param alignedCluster Cluster to be compared with all members of
     * searchedClusters
     * @param searchedClusters All members of this list will be compared with
     * alignedCluster
     * @param idString ID to be added to temporal files. This String MUST be
     * unique among all parallel instances of SingleThreadHHsearchRunner.
     */
    public SingleThreadHHsearchRunner(Cluster alignedCluster, Collection<Cluster> searchedClusters, String idString) {
        this.alignedCluster = alignedCluster;
        this.searchedClusters = searchedClusters;
        this.idString = idString;
    }

    /**
     * Runs hhsearch and returns results
     *
     * @return List of all identified matches
     * @throws Exception
     */
    @Override
    public List<HHalignHit> call() throws Exception {
        runHHsearch();
        return parseResult();
    }

    /**
     * Sets appropriate parameters and launches hhsearch
     *
     * @throws Exception
     */
    private void runHHsearch() throws Exception {
        String clusterListPath = Settings.getInstance().getHhsearchOutDirectory() + idString + ".pal";
        FileIOManager.saveClusterHHPaths(searchedClusters, clusterListPath);
        List<String> parameters = new ArrayList<>();
        parameters.add("-i");
        parameters.add(Settings.getInstance().getHhDirectory() + alignedCluster.getId() + ".hhm");
        parameters.add("-d");
        parameters.add(clusterListPath);
        parameters.add("-o");
        parameters.add(Settings.getInstance().getHhsearchOutDirectory() + idString + ".res");
        parameters.add("-cpu");
        parameters.add("1");
        if (Settings.getInstance().getHhsearchParameters() != null) {
            parameters.addAll(Settings.getInstance().getHhsearchParameters());
        }
        ExternalProcessRunner.runProcess(Settings.getInstance().getHhsearchCommand(), parameters, null, System.err, "While running hhsearch for cluster " + alignedCluster.getId() + " ", Hammock.hhSuiteEnv);
    }

    /**
     * Parses hhsearch output file.
     *
     * @return List of all hits identified by hhsearch
     * @throws IOException
     */
    private List<HHalignHit> parseResult() throws IOException {
        Map<Integer, Integer> lengthMap = null;
        Integer alignedClusterLength = null;
        if (Hammock.relativeHHScore) {
            alignedClusterLength = HHsuiteRunner.getHHLength(alignedCluster);
            lengthMap = new HashMap<>();
            for (Cluster cl : searchedClusters) {
                lengthMap.put(cl.getId(), HHsuiteRunner.getHHLength(cl));
            }
        }
        List<HHalignHit> result = new ArrayList<>();
        Map<Integer, Cluster> clusterMap = new HashMap<>();
        for (Cluster cl : searchedClusters) {
            clusterMap.put(cl.getId(), cl);
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(
                Settings.getInstance().getHhsearchOutDirectory() + idString + ".res")))) {
            String line;
            int id;
            double score;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    id = Integer.parseInt(line.substring(1));
                    line = reader.readLine();
                    String[] splitLine = line.split("\\s+"); //split on any number of whitespace characters
                    score = Double.parseDouble(splitLine[2].split("=")[1]);
                    if (Hammock.relativeHHScore) {
                        score = (score) / Math.min(alignedClusterLength, lengthMap.get(id));
                    }
                    reader.readLine();
                    String firstLine = reader.readLine();
                    reader.readLine();
                    String secondLine = reader.readLine();
                    result.add(new HHalignHit(alignedCluster, clusterMap.get(id), score, firstLine, secondLine));
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return result;
    }
}
