/* 
 * Class performs HMM-related tasks by running Hmmer as external program.
 */
package hammock;

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

/**
 *
 * @author Adam Krejci
 */
public class HmmerRunner {

    /**
     * Builds hmms for collection of clusters using parallel execution of
     * multiple instances of hmmbuild. Clusters with hasHMM == ture will be
     * ommited form hmm building.
     *
     * Sets hasHMM to true to every cluster if successful. Clusters must be
     * modifiable. Writes out files with hmms to directory specified by
     * Settings. File name is in format: [cluster_id].hmm If such file exists
     * for a cluster having hasHMM() set to false, file will be overwritten by a
     * new file without warning.
     *
     * @param clusters Clusters for which hmms will be built.
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static void buildHmms(Collection<Cluster> clusters) throws InterruptedException, ExecutionException, Exception {
        List<Callable<Void>> callables = new ArrayList<>();
        for (Cluster cluster : clusters) {
            callables.add(new SingleThreadHmmbuildRunner(cluster));
        }
        ParallelRunner.ExecuteInParallel(callables, Hammock.threadPool);
    }

    /**
     * Builds hmm for single cluster. If the cluster already hasHMM(), nothing
     * happens. Sets hasHMM to true to this cluster if successful. Cluster must
     * be modifiable. Writes out file with hmms to directory specified by
     * Settings. File name is in format: [cluster_id].hmm If such file exists
     * for a cluster having hasHMM() set to false, file will be overwritten by a
     * new file with no warning.
     *
     * @param cluster Cluster for which hmm will be built.
     * @throws IOException
     * @throws InterruptedException
     */
    public static void buildHmm(Cluster cluster) throws IOException, InterruptedException, DataException {
        SingleThreadHmmbuildRunner runner = new SingleThreadHmmbuildRunner(cluster);
        runner.call();
    }

    /**
     * Performs database search for a collection of Clusters running hmmsearch
     * as external process in multiple parallel threads. Constructs HMMs for
     * clusters if necessary. Clusters must be modifiable. Writes out database
     * sequences to file specified by Settings. Writes out search results to
     * directory specified by Settings. File name format is [cluster_id].out. If
     * such file already existed before method invocation, it will be
     * overwritten without warning. Searches aginst Parses results and returns
     * list of significant sequence hits.
     *
     * @param clusters Clusters to search the database with
     * @param sequences Seqeucnes to form the database of.
     * @return List of all hits reported by hmmsearch
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     */
    public static List<HmmsearchSequenceHit> searchWithHmms(Collection<Cluster> clusters, Collection<UniqueSequence> sequences) throws IOException, InterruptedException, ExecutionException {
        Map<String, UniqueSequence> sequenceMap = FileIOManager.saveUniqueSequencesToFastaWithoutLabels(sequences, Settings.getInstance().getFastaDatabaseFile());
        Map<Integer, Cluster> clusterMap = new HashMap<>();
        for (Cluster cluster : clusters) {
            clusterMap.put(cluster.getId(), cluster);
        }
        List<HmmsearchResult> searchResult = runParallelHmmsearch(clusters);
        List<HmmsearchSequenceHit> hits = new ArrayList<>();
        for (HmmsearchResult r : searchResult) {
            if (r != null) {
                UniqueSequence sequence = sequenceMap.get(r.getSeqId());
                double score = r.getScore();
                if (Hammock.relativeHmmScore){
                    score = (score) / Math.min(sequence.getSequence().length, r.getHmmLength());
                }
                hits.add(new HmmsearchSequenceHit(clusterMap.get(r.getHmmId()), sequence, score, r.getEvalue()));
            }
        }
        return hits;

    }

    /**
     * Executes hmmsearch in parallel (uses current temporal fasta file with
     * UniqueSequnces)
     *
     * @param clusters
     * @return
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws IOException
     */
    private static List<HmmsearchResult> runParallelHmmsearch(Collection<Cluster> clusters) throws InterruptedException, ExecutionException, IOException {
        List<HmmsearchResult> results = new ArrayList<>();
        CompletionService<List<HmmsearchResult>> resultPool = new ExecutorCompletionService<>(Hammock.threadPool);
        for (Cluster cluster : clusters) {
            resultPool.submit(new SingleThreadHmmsearchRunner(cluster));
        }
        for (int i = 0; i < clusters.size(); i++) {
            List<HmmsearchResult> res = resultPool.take().get();
            if (res != null) {
                results.addAll(res);
            }
        }
        return results;
    }
}

/**
 * Runs single instance of hmmbuild on single thread as external process.
 *
 * @author Adam Krejci
 */
class SingleThreadHmmbuildRunner implements Callable<Void> {

    private final Cluster cluster;

    /**
     * Constructor
     *
     * @param cluster Cluster for which hmm will be built.
     */
    public SingleThreadHmmbuildRunner(Cluster cluster) {
        this.cluster = cluster;
    }

    /**
     * Builds a HMM for single cluster. If the cluster does not have a MSA, MSA
     * will be build too. Sets the clusters hasMSA to true, if successful. If
     * the cluster already has HMM, nothing happens. Cluster must be modifiable.
     * Files named [cluster_id].hmm will be written in appropriate folder
     * according to Settings. If cluster already has .hmm file, but hasHMM()
     * returns false, new .hmm file will be produced and old file will be
     * overwritten without warning.
     *
     * @return
     * @throws InterruptedException
     * @throws IOException
     */
    @Override
    public Void call() throws InterruptedException, IOException, DataException {
        if (cluster.hasHMM()) {
            return null;
        }
        ClustalRunner.multipleAlignment(cluster);
        HHsuiteRunner.alnToA2M(cluster);
        List<String> parameters = new ArrayList<>();
        List<String> hmmbuildParameters = Settings.getInstance().getHmmbuildParameters();
        if (hmmbuildParameters != null) {
            parameters.addAll(hmmbuildParameters);
        }
        parameters.add(Settings.getInstance().getHmmDirectory() + cluster.getId() + ".hmm");
        parameters.add(Settings.getInstance().getMsaDirectory() + cluster.getId() + ".a2m");
        ExternalProcessRunner.runProcess(Settings.getInstance().getHmmbuildCommand(), parameters, null, System.err, "While building Hmmer hmm for cluster " + cluster.getId() + " ");
        cluster.setAsHasHMM();
        return null;
    }
}

/**
 * Runs single instance of hmmsearch on single thread as external process.
 *
 * @author Adam Krejci
 */
class SingleThreadHmmsearchRunner implements Callable<List<HmmsearchResult>> {

    private final Cluster cluster;

    public SingleThreadHmmsearchRunner(Cluster cluster) {
        this.cluster = cluster;
    }

    /**
     * Performs hmmsearch using a cluster and database. If the cluster does not
     * have hmm constructed, will construct this hmm. Database file specified in
     * Settings must exits. Writes results of database search in file
     * [cluster_id].out in folder specified by Settings. If this folder already
     * contains file witch such name, this file will be overwritten without
     * warning.
     *
     * @return List of significant sequence hits found in the database while
     * searching with this cluster.
     * @throws InterruptedException
     * @throws IOException
     * @throws DataException If database file does not exist.
     */
    @Override
    public List<HmmsearchResult> call() throws Exception {
        executeHmmsearch(cluster);
        return parseResult(cluster);
    }

    private void executeHmmsearch(Cluster cluster) throws DataException, IOException, InterruptedException {
        if (new File(Settings.getInstance().getFastaDatabaseFile()).exists() == false) {
            throw new DataException("Error. Trying to run hmmsearch with non-existing database file. "
                    + "Database file path: " + Settings.getInstance().getFastaDatabaseFile());
        }
        HmmerRunner.buildHmm(cluster);
        List<String> parameters = new ArrayList<>();
        List<String> hmmsearchParameters = Settings.getInstance().getHmmsearchParameters();
        if (hmmsearchParameters != null) {
            parameters.addAll(hmmsearchParameters);
        }
        parameters.add("--tblout");
        parameters.add(Settings.getInstance().getHmmsearchOutDirectory() + cluster.getId() + ".out");
        parameters.add(Settings.getInstance().getHmmDirectory() + cluster.getId() + ".hmm");
        parameters.add(Settings.getInstance().getFastaDatabaseFile());
        ExternalProcessRunner.runProcess(Settings.getInstance().getHmmsearchCommand(), parameters, null, System.err, "While searching sequence database with hmm " + cluster.getId() + " ");
    }

    /**
     * Parses temporal hmmsearch result file for appropriate to specified
     * cluster
     *
     * @param cluster
     * @return
     * @throws IOException
     */
    private List<HmmsearchResult> parseResult(Cluster cluster) throws IOException {
        Integer hmmLength = null;
        if (Hammock.relativeHmmScore){
            hmmLength = getHmmLength(cluster);
        }
        List<HmmsearchResult> results = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(
                Settings.getInstance().getHmmsearchOutDirectory() + cluster.getId() + ".out")))) {
            String line;
            while (((line = reader.readLine()) != null)) {
                if (!(line.startsWith("#"))) {
                    String[] splitLine = line.split("\\s+"); //split on any number of whitespace characters
                    Double score = Double.parseDouble(splitLine[5]);
                    Double evalue = Double.parseDouble(splitLine[4]);
                    HmmsearchResult lineResult = new HmmsearchResult(cluster.getId(), splitLine[0], score, hmmLength, evalue);
                    results.add(lineResult);
                }
            }

        } catch (IOException e) {
            throw new IOException(e);
        }
        if (results.size() > 0) {
            return results;
        } else {
            return null;
        }
    }

    private int getHmmLength(Cluster cl) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(
                Settings.getInstance().getHmmDirectory() + cl.getId() + ".hmm")))) {
            String line = null;
            for (int lineNum = 0; lineNum < 3; lineNum++) {
                line = reader.readLine();
            }
            String[] splitLine = line.split("\\s+"); //split on any number of whitespace characters
            return Integer.decode(splitLine[1]);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }
}

/**
 * Class represents result of one run of hmmsearch. Contains hmm and sequence id
 * to be later converted into HmmsearchSequenceHit containing the actual objects
 * @author Adam Krejci
 */
class HmmsearchResult {

    private final int hmmId;
    private final String seqId;
    private final double score;
    private final Integer hmmLength;
    private final Double evalue;

    public HmmsearchResult(int hmmId, String seqId, double score, Integer hmmLength) {
        this.hmmId = hmmId;
        this.seqId = seqId;
        this.score = score;
        this.hmmLength = hmmLength;
        this.evalue = null;
    }

    public HmmsearchResult(int hmmId, String seqId, double score, Integer hmmLength, Double evalue) {
        this.hmmId = hmmId;
        this.seqId = seqId;
        this.score = score;
        this.hmmLength = hmmLength;
        this.evalue = evalue;
    }

    public int getHmmId() {
        return hmmId;
    }

    public String getSeqId() {
        return seqId;
    }

    public double getScore() {
        return score;
    }

    public Integer getHmmLength() {
        return hmmLength;
    }
    
    public Double getEvalue() {
        return evalue;
    }
}
