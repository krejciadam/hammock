/*
 * Class handles most of I/O work with files
 */
package cz.krejciadam.hammock;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Adam Krejci
 */
public class FileIOManager {

    /**
     * Loads a scoring matrix. Matrix uses standard amino acid ordering (defined
     * in UniqueSequence). Input file must contain scoring matrix in BioJava
     * format.
     *
     * @param matrixFilePath path to file with scoring matrix in BioJava format
     * @return two dimensional array representing scoring matrix.
     * @throws IOException
     * @throws HammockException
     */
    public static int[][] loadScoringMatrix(String matrixFilePath) throws IOException, HammockException {
        int[][] scoringMatrix = new int[24][24];

        try (BufferedReader reader = new BufferedReader(new FileReader(new File(matrixFilePath)))) {
            String line;
            int lineCounter = 0;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("\\s")) {
                    if (!(line.replaceAll("\\s+", "").equals(new String(UniqueSequence.getAminoAcidsAndSpecials())))) {
                        throw new FileFormatException("Error in scoring matrix file: " + matrixFilePath + ". Scoring"
                                + "matrix should have exactly this order of rows/columns: " + new String(UniqueSequence.getAminoAcidsAndSpecials()));
                    }
                }
                if (!(line.startsWith("#") || line.startsWith(" ") || line.startsWith("\t"))) { //skip comments and AA header line
                    String[] splitLine = line.split("\\s+"); //split on any number of whitespace characters
                    if (splitLine.length != 25) {
                        throw new FileFormatException("Error in scoring matrix file: " + matrixFilePath + ". Scoring matrix "
                                + "should always have 24 columns (plus 1 column describing AAs).");
                    }
                    for (int i = 1; i < splitLine.length; i++) { //from 1 - first is AA name.
                        scoringMatrix[lineCounter][i - 1] = Integer.parseInt(splitLine[i]);
                    }
                    lineCounter++;
                    if (lineCounter > 24) {
                        throw new FileFormatException("Error in scoring matrix file: " + matrixFilePath + ". Scoring matrix "
                                + "should always have 24 rows (plus 1 column describing AAs).");
                    }
                }

            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new FileFormatException("Error in scoring matrix file: " + matrixFilePath + ". Scoring matrix "
                    + "should always have 24 rows (plus 1 column describing AAs).");
        }
        return scoringMatrix;
    }

    /**
     * Loads frequency matrix of background frequencies for KLD calculation
     *
     * @param path A path to the matrix
     * @return The frequency matrix
     * @throws IOException
     */
    public static Map<Character, Map<Character, Double>> loadFrequencyMatrix(String path) throws IOException {
        Map<Character, Map<Character, Double>> result = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(path)))) {
            String line;
            List<Character> aaList = new ArrayList<>();
            int lineIndex = 0;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                if (line.matches("^[A-Z]+.*")) { //start line
                    for (String s : line.split("\t")) {
                        aaList.add(s.charAt(0));
                    }
                    continue;
                }
                String[] numberSplit = line.split("\t");
                Map<Character, Double> aaMap = new HashMap<>();
                for (int i = 0; i < numberSplit.length; i++) {
                    aaMap.put(aaList.get(i), Double.parseDouble(numberSplit[i]));
                }
                result.put(aaList.get(lineIndex), aaMap);
                lineIndex++;
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return (result);
    }

    /**
     * Loads a matrix of empirical probabilities
     *
     * @param path
     * @return
     * @throws IOException
     */
    public static Map<Double, Double> loadEmpiricalProbabs(String path) throws IOException {
        Map<Double, Double> result = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(path)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] splitLine = line.split("\t");
                Double score = Double.parseDouble(splitLine[0]);
                Double probab = Double.parseDouble(splitLine[1]);
                result.put(score, probab);
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return (result);
    }

    /**
     *     /**
     * Loads a fasta file. Retains input order (of first occurrence of a
     * seuquence) Each seuqence must have a header line starting with ">".
     * Header may contain sequence label in format: ">any_text|count|label". To
     * each sequence not having label in the header, default label "no_label" is
     * assigned. Sequences in fasta file need not to be unique. Multi-line fasta
     * files are not supported.
     *
     * @param fileName Path to file to be loaded
     * @return List of UniqueSequence objects parsed in the input order (of
     * first occurrence)
     * @throws IOException
     * @throws FileNotFoundException
     * @throws FileFormatException
     */
    public static List<UniqueSequence> loadUniqueSequencesFromFasta(String fileName) throws IOException, FileNotFoundException, FileFormatException {
        Map<String, Map<String, Integer>> sequenceMap = new LinkedHashMap<>();
        String line;
        String sequence = "";
        String label = null;
        Integer count = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (sequence.length() > 0) { //add the previous sequence
                        Map<String, Integer> labelsMap = updateLabelsMap(sequenceMap.get(sequence), label, count);
                        sequenceMap.put(sequence, labelsMap);
                        sequence = "";
                    }
                    String[] splitLine = line.trim().substring(1).split("\\|");
                    if (splitLine.length >= 2) {
                        count = Integer.decode(splitLine[1].trim());
                        if (count < 1) {
                            throw new FileFormatException("Error while loading input file. Fasta header defines sequence count lower than 1.");
                        }
                    } else {
                        count = 1;
                    }
                    if (splitLine.length >= 3) {
                        label = splitLine[2];
                    } else {
                        label = "no_label";
                    }
                } else {
                    if (label == null || count == null) {
                        throw new FileFormatException("Error. Incorrect fasta format. Maybe header or sequence line missing?");
                    }
                    sequence = sequence.concat(line.trim());
                }
            }//Add the last sequence:
            Map<String, Integer> labelsMap = updateLabelsMap(sequenceMap.get(sequence), label, count);
            sequenceMap.put(sequence, labelsMap);
        }
        List<UniqueSequence> result = new ArrayList<>();
        for(Map.Entry<String, Map<String, Integer>> entry : sequenceMap.entrySet()){
            result.add(new UniqueSequence(entry.getKey(), entry.getValue()));
        }
        return (result);
    }

    private static Map<String, Integer> updateLabelsMap(Map<String, Integer> labelsMap, String label, int count) {
        if (labelsMap == null) {
            labelsMap = new HashMap<>();
            labelsMap.put(label, count);
        } else {
            Integer oldCount = labelsMap.get(label);
            if (oldCount == null) {
                oldCount = 0;
            }
            labelsMap.put(label, oldCount + count);
        }
        return (labelsMap);
    }

    /**
     * Loads UniqueSequence from a .csv table. Table must contain header
     * specifying full list of labels.
     *
     * @param fileName
     * @return
     * @throws IOException
     * @throws cz.krejciadam.hammock.FileFormatException if wrong amino acid letter is used
     */
    public static List<UniqueSequence> loadUniqueSequencesFromTable(String fileName) throws IOException, FileFormatException {
        List<UniqueSequence> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line = reader.readLine(); //header line
            String[] headerLine = line.split(Hammock.CSV_SEPARATOR);
            List<String> labels = new ArrayList<>();
            for (int i = 1; i < headerLine.length; i++) {
                labels.add(headerLine[i]);
            }
            while ((line = reader.readLine()) != null) {
                List<String> splitLine = Arrays.asList(line.split(Hammock.CSV_SEPARATOR));
                result.add(lineToUniqueSequence(splitLine.get(0), splitLine.subList(1, splitLine.size()), labels));
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return result;
    }

    private static UniqueSequence lineToUniqueSequence(String sequence, List<String> labelsLine, List<String> labels) throws FileFormatException {
        Map<String, Integer> labelsMap = new HashMap<>();
        for (int i = 0; i < labelsLine.size(); i++) {
            int value = Integer.decode(labelsLine.get(i));
            if (value != 0) {
                labelsMap.put(labels.get(i), value);
            }
        }
        return (new UniqueSequence(sequence, labelsMap));
    }

    /**
     * For specified clusters, loads their alignments from a file, saves them as
     * temporal files and marks the clusters as having MSAs. If some clusters
     * are missing or different in size in the file, a FileFormatException is
     * thrown.
     *
     * @param clusters Clusters whose alignments are to be loaded
     * @param fileName Path to the file with saved clusters
     * @return
     * @throws FileFormatException
     * @throws IOException
     */
    public static List<Cluster> loadClusterAlignmentsFromFile(List<Cluster> clusters, String fileName) throws FileFormatException, IOException {
        List<Cluster> result = new ArrayList<>();
        Set<Integer> ids = new HashSet<>();
        for (Cluster cl : clusters) {
            ids.add(cl.getId());
        }
        List<Cluster> loadedClusters = loadClusterDetailsFromCsv(fileName, ids, false);
        Map<Integer, Cluster> loadedClustersMap = new HashMap<>();
        for (Cluster cl : loadedClusters) {
            loadedClustersMap.put(cl.getId(), cl);
        }
        for (Cluster cl : clusters) {
            if (cl.size() != loadedClustersMap.get(cl.getId()).size()) {
                throw new IOException("Cluster " + cl.getId() + " has different size than the number of alignment lines in file " + fileName);
            }
            result.add(loadedClustersMap.get(cl.getId()));
        }
        return result;
    }

    /**
     * Loads clusters from a .csv file. The file may or may not contain
     * alignment column. Clusters will be loaded with or without alignments,
     * accordingly.
     *
     * @param fileName A path to the clusters' csv file.
     * @param loadAlignments Should alignments be loaded from the file?
     * @return
     * @throws FileFormatException
     * @throws IOException
     */
    public static List<Cluster> loadClustersFromCsv(String fileName, boolean loadAlignments) throws FileFormatException, IOException {
        return (loadClusterDetailsFromCsv(fileName, null, loadAlignments));
    }

    private static List<Cluster> loadClusterDetailsFromCsv(String fileName, Set<Integer> alignmentsToGet, boolean loadAllAlignments) throws FileFormatException, IOException {
        List<Cluster> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            List<String> header = new ArrayList<>(Arrays.asList(reader.readLine().split(Hammock.CSV_SEPARATOR)));
            int alignmentIndex = header.indexOf("alignment");
            if (alignmentIndex != -1) { //to correct the next index for removing this
                header.remove(alignmentIndex);
            }
            int sumIndex = header.indexOf("sum");
            if (sumIndex != -1) {
                header.remove(sumIndex);
            }
            List<String> labels = header.subList(2, header.size()); //0: cluster_id, 1: sequence
            String line;
            List<UniqueSequence> sequences = null;
            StringBuilder alignments = new StringBuilder();
            int clusterId = -1;
            int seqIndex = 0;
            while ((line = reader.readLine()) != null) {
                List<String> splitLine = new ArrayList<>(Arrays.asList(line.split(Hammock.CSV_SEPARATOR)));
                int id = Integer.decode(splitLine.get(0));
                if (id != clusterId) {
                    if (clusterId != -1) {  //it is not the first cluster => add the cluster
                        Cluster cl = new Cluster(sequences, clusterId);
                        if ((alignments.length() >= 1) && ((loadAllAlignments) || ((alignmentsToGet != null) && (alignmentsToGet.contains(clusterId))))) {
                            int alignmentLinesCount = (StringUtils.countMatches(alignments, "\n") / 2);
                            if (alignmentLinesCount == cl.getUniqueSize()) { //alignments with some aligned sequences 'NA' will not be used
                                FileIOManager.saveStringToFile(alignments.toString(), Settings.getInstance().getMsaDirectory() + clusterId + ".aln");
                                cl.setAsHasMSA();
                            }
                        }
                        result.add(cl);
                    }
                    clusterId = id;
                    sequences = new ArrayList<>();
                    alignments = new StringBuilder();
                    seqIndex = 0;
                }
                if (alignmentIndex != -1) {
                    String alignedSeq = splitLine.get(alignmentIndex);
                    if (!alignedSeq.equals("NA")) {
                        alignments.append(">").append(clusterId).append("_").append(seqIndex).append("\n").append(alignedSeq).append("\n");
                    }
                    splitLine.remove(alignmentIndex);
                    seqIndex++;
                }
                if (sumIndex != -1) {
                    splitLine.remove(sumIndex);
                }
                sequences.add(lineToUniqueSequence(splitLine.get(1), splitLine.subList(2, splitLine.size()), labels));
            }
            Cluster cl = new Cluster(sequences, clusterId);
            if (alignments.length() >= 1) {
                FileIOManager.saveStringToFile(alignments.toString(), Settings.getInstance().getMsaDirectory() + clusterId + ".aln");
                cl.setAsHasMSA();
            }
            result.add(cl);
        } catch (NumberFormatException | NullPointerException e) {
            throw new FileFormatException("Error in cluster file: "
                    + fileName + " - wrong format. Original message: ", e);
        }
        return (result);
    }

    /**
     * @see FileIOManager#saveClusterSequencesToCsvOrdered(Collection<Cluster>,
     * String, List<String>, Collection<UniqueSequence>, Map<Cluster, String>)
     */
    public static void saveClusterSequencesToCsvOrdered(Collection<Cluster> clusters, String filePath, List<String> labels, Collection<UniqueSequence> orderedSequences) throws IOException {
        Map<Cluster, List<String>> alignmentsMap = getAlignmentsMap(clusters);
        writeClusterSequencesToCsv(orderedSequences, clusters, filePath, labels, alignmentsMap);
    }

    /**
     * Saves a .csv file containing one line per sequence for every cluster on
     * input, retains the order defined by orderedSequences
     *
     * @param clusters Clusters to be saved
     * @param filePath Path to resulting file
     * @param labels List of labels. Only these labels and their appropriate
     * counts will be saved in resulting file.
     * @param orderedSequences Defines the order of the sequences in the file
     * @param alignmentsMap An existing alignmentsMap, which will be re-formated
     * for writing purposes
     * @throws IOException
     */
    public static void saveClusterSequencesToCsvOrdered(Collection<Cluster> clusters, String filePath, List<String> labels, Collection<UniqueSequence> orderedSequences, Map<Cluster, String> alignmentsMap) throws IOException {
        Map<Cluster, List<String>> newAlignmentsMap = reformatAlignmentsMap(clusters, alignmentsMap);
        writeClusterSequencesToCsv(orderedSequences, clusters, filePath, labels, newAlignmentsMap);
    }

    /**
     * @see saveClusterSequencesToCsv(Collection<Cluster>, String, List<String>,
     * Map<Cluster, String>)
     */
    public static void saveClusterSequencesToCsv(Collection<Cluster> clusters, String filePath, List<String> labels) throws IOException {
        List<Cluster> clusterSortedList = new ArrayList<>(clusters);
        Collections.sort(clusterSortedList, Collections.reverseOrder());
        List<UniqueSequence> sortedSequences = getSortedSequences(clusterSortedList);
        Map<Cluster, List<String>> alignmentsMap = getAlignmentsMap(clusters);
        writeClusterSequencesToCsv(sortedSequences, clusters, filePath, labels, alignmentsMap);
    }

    /**
     * Saves a csv file containing one line per sequence for every cluster on
     * input. The resulting file is ordered from the largest cluster to the
     * smallest, within a cluster, sequences are ordered by size and
     * alphabetically
     *
     * @param clusters Clusters to be saved
     * @param filePath Path to resulting file
     * @param labels List of labels. Only these labels and their appropriate
     * counts will be saved in resulting file.
     * @param alignmentsMap An existing alignment map which will be re-formated
     * for writing purposes
     * @throws IOException
     */
    public static void saveClusterSequencesToCsv(Collection<Cluster> clusters, String filePath, List<String> labels, Map<Cluster, String> alignmentsMap) throws IOException {
        List<Cluster> clusterSortedList = new ArrayList<>(clusters);
        Collections.sort(clusterSortedList, Collections.reverseOrder());
        List<UniqueSequence> sortedSequences = getSortedSequences(clusterSortedList);
        Map<Cluster, List<String>> newAlignmentsMap = reformatAlignmentsMap(clusters, alignmentsMap);
        writeClusterSequencesToCsv(sortedSequences, clusters, filePath, labels, newAlignmentsMap);
    }

    private static double empiricalProbab(double score, double minScore, double maxScore, Map<Double, Double> empiricalProbabs) {
        if (score < minScore) {
            return 1;
        }
        if (score > maxScore) {
            return 0;
        }
        return (empiricalProbabs.get((double) Math.round(score * 10) / 10)); //one decimal place rounding
    }

    /**
     * Saves the results of multiple HMMSearch runs to a csv file.
     *
     * @param hits The HMMSearch hits to be saved
     * @param filePath Path to save the results into
     * @param empiricalProbabsFile File containing empirical probabilities for
     * this search
     * @param clusterCount How many clusters are in the result
     * @param sequenceCount How many sequences are in the result
     * @throws IOException
     */
    public static void saveHmmsearchHitsToCsv(Collection<HmmsearchSequenceHit> hits, String filePath, String empiricalProbabsFile, int clusterCount, int sequenceCount) throws IOException {
        Map<Double, Double> empiricalProbabs = null;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        if (empiricalProbabsFile != null) {
            empiricalProbabs = FileIOManager.loadEmpiricalProbabs(empiricalProbabsFile);
            for (Double score : empiricalProbabs.keySet()) {
                if (score > max) {
                    max = score;
                }
                if (score < min) {
                    min = score;
                }
            }
        }
        List<HmmsearchSequenceHit> hitList = new ArrayList<>(hits);
        Collections.sort(hitList, Collections.reverseOrder());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("cluster_id" + Hammock.CSV_SEPARATOR
                    + "main_sequence" + Hammock.CSV_SEPARATOR
                    + "found_sequence" + Hammock.CSV_SEPARATOR
                    + "score");
            if (empiricalProbabs != null) {
                writer.write(Hammock.CSV_SEPARATOR + "evalue_empirical");
            }
            writer.write("\n");
            for (HmmsearchSequenceHit hit : hitList) {
                writer.write(hit.getCluster().getId() + Hammock.CSV_SEPARATOR
                        + hit.getCluster().getSequences().get(0).getSequenceString() + Hammock.CSV_SEPARATOR
                        + hit.getSequence().getSequenceString() + Hammock.CSV_SEPARATOR
                        + hit.getScore());
                if (empiricalProbabs != null) {
                    writer.write(Hammock.CSV_SEPARATOR + empiricalProbab(hit.getScore(), min, max, empiricalProbabs) * clusterCount * sequenceCount);
                }
                writer.write("\n");
            }
        }
    }

    /**
     * Saves the results of multiple HHAlign runs to a csv file in the form of
     * an all vs. all matrix.
     *
     * @param hits The HHAlign hits to be saved
     * @param clusters1 One set of compared clusters (matrix rows)
     * @param clusters2 The other set of compared clusters (matrix columns)
     * @param filePath A path to save the file to.
     * @throws IOException
     */
    public static void saveHHAlignHitsToCsv(Collection<HHalignHit> hits, Collection<Cluster> clusters1, Collection<Cluster> clusters2, String filePath) throws IOException {
        Map<Integer, Map<Integer, Double>> scoreMatrix = new HashMap<>();
        for (HHalignHit hit : hits) {
            int lineId = hit.getSearchedCluster().getId();
            int columnId = hit.getFoundCluster().getId();
            Map<Integer, Double> matrixLine;
            matrixLine = scoreMatrix.get(lineId);
            if (matrixLine == null) {
                matrixLine = new HashMap<>();
            }
            matrixLine.put(columnId, hit.getScore());
            scoreMatrix.put(lineId, matrixLine);
        }
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (Cluster cl2 : clusters2) {
                writer.write(Hammock.CSV_SEPARATOR + cl2.getId());
            }
            for (Cluster cl1 : clusters1) {
                writer.newLine();
                writer.write("" + cl1.getId());
                for (Cluster cl2 : clusters2) {
                    Double score = scoreMatrix.get(cl1.getId()).get(cl2.getId());
                    if (score == null) {
                        score = -1.0;
                    }
                    writer.write(Hammock.CSV_SEPARATOR + score);
                }
            }
        }
    }

    private static List<UniqueSequence> getSortedSequences(List<Cluster> sortedClusters) {
        List<UniqueSequence> sortedSequences = new ArrayList<>();
        for (Cluster cl : sortedClusters) {
            List<UniqueSequence> sequences = cl.getSequences();
            Collections.sort(sequences, Collections.reverseOrder(new UniqueSequenceSizeAlphabeticComparator()));
            sortedSequences.addAll(sequences);
        }
        return (sortedSequences);
    }

    private static Map<Cluster, List<String>> reformatAlignmentsMap(Collection<Cluster> clusters, Map<Cluster, String> alignmentsMap) {
        Map<Cluster, List<String>> newAlignmentsMap = new HashMap<>();
        for (Cluster cl : clusters) {
            if (alignmentsMap.containsKey(cl)) {
                List<String> alignmentLines = new ArrayList<>();
                String[] split = alignmentsMap.get(cl).split("\n");
                for (int i = 1; i < split.length; i += 2) { //odd lines only
                    alignmentLines.add(split[i]);
                }
                newAlignmentsMap.put(cl, alignmentLines);
            }
        }
        return (newAlignmentsMap);
    }

    private static Map<Cluster, List<String>> getAlignmentsMap(Collection<Cluster> clusters) throws IOException {
        Map<Cluster, List<String>> alignmentsMap = new HashMap<>();
        for (Cluster cl : clusters) {
            if (cl.hasMSA()) {
                alignmentsMap.put(cl, FileIOManager.getAlignmentLines(cl));
            }
        }
        return (alignmentsMap);
    }

    /**
     * Saves a collection of unique sequences into a file
     *
     * @param sequences Sequences to be saved
     * @param filePath File to write the sequences into
     * @param labels A list of sequence labels. Defines the column order
     * @throws IOException
     */
    public static void saveUniqueSequencesToCsv(Collection<UniqueSequence> sequences, String filePath, List<String> labels) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("sequence");
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + label);
            }
            writer.newLine();
            for (UniqueSequence seq : sequences) {
                writer.write(seq.getSequenceString());
                for (String label : labels) {
                    Integer count = seq.getLabelsMap().get(label);
                    if (count == null) {
                        count = 0;
                    }
                    writer.write(Hammock.CSV_SEPARATOR + count);
                }
                writer.newLine();
            }
        }
    }

    private static void writeClusterSequencesToCsv(Collection<UniqueSequence> sequences, Collection<Cluster> clusters, String filePath, List<String> labels, Map<Cluster, List<String>> clusterAlignments) throws IOException {
        Map<String, String> msaMap = new HashMap<>();
        Map<String, Cluster> sequenceClusterMap = new HashMap<>();
        for (Cluster cl : clusters) {
            for (UniqueSequence seq : cl.getSequences()) {
                sequenceClusterMap.put(seq.getSequenceString(), cl);
            }
            if (clusterAlignments.containsKey(cl)) {
                for (String msaLine : clusterAlignments.get(cl)) {
                    msaMap.put(msaLine.replace("-", ""), msaLine);
                }
            }
        }
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("cluster_id" + Hammock.CSV_SEPARATOR + "sequence" + Hammock.CSV_SEPARATOR + "alignment" + Hammock.CSV_SEPARATOR + "sum");
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + label);
            }
            writer.newLine();
            for (UniqueSequence seq : sequences) {
                Cluster cluster = sequenceClusterMap.get(seq.getSequenceString());
                if (cluster != null) {
                    writer.write(cluster.getId() + Hammock.CSV_SEPARATOR + seq.getSequenceString() + Hammock.CSV_SEPARATOR);
                    if (msaMap.containsKey(seq.getSequenceString())) {
                        writer.write(msaMap.get(seq.getSequenceString()) + Hammock.CSV_SEPARATOR);
                    } else {
                        writer.write("NA" + Hammock.CSV_SEPARATOR);
                    }
                } else {
                    writer.write("NA" + Hammock.CSV_SEPARATOR + seq.getSequenceString() + Hammock.CSV_SEPARATOR + "NA" + Hammock.CSV_SEPARATOR);
                }
                writer.write("" + seq.size());
                for (String label : labels) {
                    Integer count = seq.getLabelsMap().get(label);
                    if (count == null) {
                        count = 0;
                    }
                    writer.write(Hammock.CSV_SEPARATOR + count);
                }
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves cluster overview to .csv file. Lines are sorted by cluster size
     * from the largest.
     *
     * @param clusters Clusters to be saved
     * @param filePath Path to a file write the clusters into
     * @param labels A list of labels to use. Defines column order
     * @throws IOException
     */
    public static void SaveClustersToCsv(Collection<Cluster> clusters, String filePath, List<String> labels) throws IOException {
        List<Cluster> sortedList = new ArrayList<>(clusters);
        Collections.sort(sortedList, Collections.reverseOrder());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("cluster_id" + Hammock.CSV_SEPARATOR + "main_sequence" + Hammock.CSV_SEPARATOR + "sum");
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + label);
            }
            writer.newLine();
            for (Cluster cl : sortedList) {
                List<UniqueSequence> sequences = cl.getSequences();
                Collections.sort(sequences, Collections.reverseOrder());
                writer.write(cl.getId() + Hammock.CSV_SEPARATOR + sequences.get(0).getSequenceString() + Hammock.CSV_SEPARATOR + cl.size());
                Map<String, Integer> clusterCountMap = getClusterLabelsMap(cl);
                for (String label : labels) {
                    Integer count = clusterCountMap.get(label);
                    if (count == null) {
                        count = 0;
                    }
                    writer.write(Hammock.CSV_SEPARATOR + count);
                }
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException(e);
        }

    }

    /**
     * Returns labels count map of an entire cluster (sum of all sequences)
     *
     * @param cl
     * @return
     */
    private static Map<String, Integer> getClusterLabelsMap(Cluster cl) {
        Map<String, Integer> res = new HashMap<>();
        for (UniqueSequence seq : cl.getSequences()) {
            for (Map.Entry<String, Integer> entry : seq.getLabelsMap().entrySet()) {
                Integer current = res.get(entry.getKey());
                if (current == null) {
                    current = 0;
                }
                current += entry.getValue();
//                current += 1;
                res.put(entry.getKey(), current);
            }
        }
        return res;
    }

    /**
     * Saves a file containing input statistics - information on sequences and
     * their counts
     *
     * @param sequences
     * @param labels Labels to be part of the file. Defines column order
     * @param filePath
     * @throws IOException
     */
    public static void saveInputStatistics(Collection<UniqueSequence> sequences, List<String> labels, String filePath) throws IOException {
        Map<String, Integer> totalMap = getTotalLabelCounts(sequences, labels);
        Map<String, Integer> uniqueMap = getUniqueLabelCounts(sequences, labels);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + label);
            }
            writer.newLine();
            writer.write("total_count");
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + totalMap.get(label));
            }
            writer.newLine();
            writer.write("unique_count");
            for (String label : labels) {
                writer.write(Hammock.CSV_SEPARATOR + uniqueMap.get(label));
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Loads the alignments from a MSA file.
     *
     * @param path Path to the MSA file
     * @return A list of strings, each member is one MSA string. Input order is
     * retained.
     * @throws IOException
     */
    public static List<String> getAlignmentLines(String path) throws IOException {
        List<String> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(path)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = reader.readLine();               //only even lines
                result.add(line.trim());
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return result;
    }

    /**
     * Loads temporal MSA file for a cluster and returns only lines containing
     * alignment in form of a (multiline) string
     *
     * @param cl
     * @return
     * @throws IOException
     */
    public static List<String> getAlignmentLines(Cluster cl) throws IOException {
        List<String> list = new ArrayList<>(Arrays.asList(getAllAlignmentLines(cl)));
        List<String> res = new ArrayList<>();
        for (int i = 1; i < list.size(); i += 2) {
            res.add(list.get(i));
        }
        return (res);
    }

    private static String[] getAllAlignmentLines(Cluster cl) throws IOException {
        if (cl.getUniqueSize() > 1) {
            return (fileAsString(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln").split("\n"));
        } else {
            return (cl.getFastaString().split("\n"));
        }
    }

    /**
     * Saves one or more clusters to one fasta File.
     *
     * @param clusters clusters to be saved
     * @param filePath path to output .fasta file
     * @throws IOException
     */
    public static void saveClustersToFasta(Collection<Cluster> clusters, String filePath) throws IOException {
        String toWrite = "";
        for (Cluster cluster : clusters) {
            toWrite = toWrite.concat(cluster.getFastaString());
        }
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(toWrite);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves one Cluster to a fasta file. No linebreaks are addeded. Cluster's
     * sequences are saved in reversed natural order.
     *
     * @param cluster Cluster object to be saved
     * @param filePath Path to a file Cluster objects will be saved to.
     * @throws IOException
     */
    public static void saveClusterToFasta(Cluster cluster, String filePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(cluster.getFastaString());
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Returns the first aligned sequence (including gaps) from a MSA in .aln
     * format (aligned fasta)
     *
     * @param cl cluster to get the first aligned sequence for
     * @return the first aligned sequence of a .aln MSA file
     * @throws DataException
     * @throws IOException
     */
    public static String getFirstAlnSeq(Cluster cl) throws DataException, IOException {
        if (!cl.hasMSA()) {
            throw new DataException("Error, Can't return alignment line for a cluster without alignment.");
        }
        return getNthLine(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln", 1);
    }

    /**
     * Returns the first line of a cluster's alignment in the A2m format.
     *
     * @param cl
     * @return
     * @throws DataException
     * @throws IOException
     */
    public static String getFirstA2mSeq(Cluster cl) throws DataException, IOException {
        if (!cl.hasMSA()) {
            throw new DataException("Error, Can't return alignment line for a cluster without alignment.");
        }
        return getNthLine(Settings.getInstance().getMsaDirectory() + cl.getId() + ".a2m", 1);
    }

    private static String getNthLine(String file, int skip) throws IOException {
        String res = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            for (int i = 0; i < skip; i++) {
                reader.readLine();
            }
            res = reader.readLine();
        } catch (IOException e) {
            throw new IOException(e);
        }
        return res;
    }

    /**
     * Merges MSA files for two clusters, adds any gap columns specified
     *
     * @param cl1 first cluster to be merged
     * @param cl2 second cluster to be merged
     * @param gaps1 gaps to be added into first cluster
     * @param gaps2 gaps to be added into second cluster
     * @param newId id of resulting cluster
     * @throws IOException
     */
    public static void mergeAlignedClusters(Cluster cl1, Cluster cl2, List<Integer> gaps1, List<Integer> gaps2, int newId) throws IOException {
        String cl1String = insertGapsIntoAlignment(cl1, gaps1, newId, 1);
        String cl2String = insertGapsIntoAlignment(cl2, gaps2, newId, cl1.getUniqueSize() + 1);
        cl1String = cl1String.concat(cl2String);
        saveStringToFile(cl1String, Settings.getInstance().getMsaDirectory() + newId + ".aln");
    }

    /**
     * Reformats a MSA file. Adds gaps on positions specified by gapPositions
     * parameter
     *
     * @param id cluster to which .aln MSA file to add gaps
     * @param gapPositions list of position of new gaps
     * @param newId id of resulting .aln file
     * @param startingSequenceId first sequence will have this id, id will
     * increade throughout the file
     * @return String reprecenting whole MSA in .aln format
     * @throws IOException
     */
    private static String insertGapsIntoAlignment(Cluster cl, List<Integer> gapPositions, int newId, int startingSequenceId) throws IOException {
        StringBuilder res = new StringBuilder();
        int currentId = startingSequenceId;
        String[] alnLines = getAllAlignmentLines(cl);
        for (String line : alnLines) {
            if (line.startsWith(">")) {
                res.append(">").append(newId).append("_").append(currentId).append("\n");
                currentId++;
            } else {
                StringBuilder lineBuilder = new StringBuilder(line);
                for (int gapPosition : gapPositions) {
                    lineBuilder.insert(gapPosition, '-');
                }
                res.append(lineBuilder).append("\n");
            }
        }
        return res.toString();
    }

    /**
     * Creates an alignment out of a greedy incremental clustering result
     *
     * @param pivot The central sequence of the greedy cluster
     * @param results Sequences aligned to the pivot
     * @param clusterId The cluster id to be part of the result
     * @return A string containing the alignment in the form of a .aln (aligned
     * fasta) file.
     * @throws IOException
     */
    public static String makeShiftedClusterAlignment(UniqueSequence pivot, List<AligningScorerResult> results, int clusterId) throws IOException {
        int maxLengthPlusShift = pivot.getSequence().length; //if everything is shifted to left or short, pivot is the longest
        int minShift = 0;
        for (AligningScorerResult result : results) {
            if (result.getShift() < minShift) {
                minShift = result.getShift();
            }
            if ((result.getSequence().getSequence().length + result.getShift()) > maxLengthPlusShift) {
                maxLengthPlusShift = result.getSequence().getSequence().length + result.getShift();
            }
        }
        StringBuilder alignment = new StringBuilder();
        int index = 0;
        alignment.append(">").append(clusterId).append("_").append(index).append("\n");
        for (int i = 0; i < Math.abs(minShift); i++) {
            alignment.append('-');
        }
        alignment.append(pivot.getSequenceString());
        for (int i = 0; i < (maxLengthPlusShift - pivot.getSequence().length); i++) {
            alignment.append('-');
        }
        alignment.append("\n");
        for (AligningScorerResult result : results) {
            index++;
            alignment.append(">").append(clusterId).append("_").append(index).append("\n");
            int actualGapsLeft;
            int actualGapsRight;
            if (result.getShift() <= 0) {
                actualGapsLeft = Math.abs(minShift) - Math.abs(result.getShift());
                actualGapsRight = maxLengthPlusShift - result.getSequence().getSequence().length + (Math.abs(minShift) - actualGapsLeft);
            } else {
                actualGapsRight = maxLengthPlusShift - (result.getSequence().getSequence().length + result.getShift());
                actualGapsLeft = Math.abs(minShift) + result.getShift();
            }
            for (int i = 0; i < actualGapsLeft; i++) {
                alignment.append('-');
            }
            alignment.append(result.getSequence().getSequenceString());
            for (int i = 0; i < actualGapsRight; i++) {
                alignment.append('-');
            }
            if (index <= results.size()) {
                alignment.append('\n');
            }
            if (pivot.getSequenceString().equals(result.getSequence().getSequenceString())) {
                throw new IOException("Error. Duplicate pivot.");
            }
        }
        return (alignment.toString());
    }

    /**
     * For every cluster, saves the path to its appropriate .hmm file in
     * temporal files
     *
     * @param clusters
     * @param filePath
     * @throws IOException
     */
    public static void saveClusterHHPaths(Collection<Cluster> clusters, String filePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (Cluster cl : clusters) {
                writer.write(Settings.getInstance().getHhDirectory() + Settings.getInstance().getSeparatorChar() + cl.getId() + ".hhm");
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Writes a collection of sequences intoto a fasta file. Assigns a unique
     * (unique only within this file) id to each sequences, this id will be
     * written to fasta file. Returns a map where key is assigned unique id and
     * value is UniqueSequence object with this unique id assigned. Does not
     * save information on unique sequence's labels and label counts.
     *
     * Each Cluster's sequences are saved in reversed natural order.
     *
     * @param sequences UniqueSequence objects to be written to fasta file
     * @param filePath path to a file UniqeSequence objects will be saved to
     * @return A map mapping new unique id used in output fasta file to each
     * UniqueSequence.
     * @throws IOException
     */
    public static Map<String, UniqueSequence> saveUniqueSequencesToFastaWithoutLabels(Collection<UniqueSequence> sequences, String filePath) throws IOException {
        Map<String, UniqueSequence> idMap = new HashMap<>();
        int id = 0;
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (UniqueSequence sequence : sequences) {
                id++;
                writer.write(">" + id);
                writer.newLine();
                writer.write(sequence.getSequenceString());
                writer.newLine();
                idMap.put("" + id, sequence);
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return idMap;
    }

    /**
     * Writes a fasta file containing the unique sequences on input. Label and
     * count information will be retained in the fasta headers. Input order is
     * retained.
     *
     * @param sequences Sequences to save
     * @param filePath Path to a file to save the sequences into.
     * @throws IOException
     */
    public static void saveUniqueSequencesToFasta(Collection<UniqueSequence> sequences, String filePath) throws IOException {
        int id = 0;
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (UniqueSequence seq : sequences) {
                for (Map.Entry<String, Integer> entry : seq.getLabelsMap().entrySet()) {
                    writer.write(">" + id + "|" + entry.getValue() + "|" + entry.getKey());
                    writer.newLine();
                    writer.write(seq.getSequenceString());
                    writer.newLine();
                    id++;
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves a single UniqueSequence object into a fasta file (containing then
     * just one record)
     *
     * @param sequence A UniqueSequence object to be saved
     * @param filePath path to file to write to.
     * @param id
     * @throws IOException
     */
    public static void saveUniqueSequenceToFastaWithoutLabels(UniqueSequence sequence, String filePath, String id) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(">" + id);
            writer.newLine();
            writer.write(sequence.getSequenceString());
            writer.newLine();
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Renames all records in a fasta file so that seuqnce headers are in format
     * >newId_n where n is sequence index (from top to bottom in file)
     *
     * @param filePath path to fasta file to be changed
     * @param newId header prefix
     * @throws IOException
     */
    public static void renameFasta(String filePath, String newId) throws IOException {
        String toWrite = "";
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            int i = 1;
            while ((line = reader.readLine()) != null) { //odd line
                toWrite = toWrite.concat(">" + newId + "_" + i + "\n");
                i++;
                toWrite = toWrite.concat(reader.readLine() + "\n"); //even line
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(toWrite);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Deletes all files in a specified folder
     *
     * @param folderName Path to folder in which all files will be deleted
     */
    public static void deleteFolderContents(String folderName) {
        File folder = new File(folderName);
        File[] files = folder.listFiles();
        if (files != null) { //some JVMs return null for empty dirs
            for (File f : files) {
                f.delete();
            }
        }
    }

    /**
     * Returns absolute paths to all files in input folder
     *
     * @param folderName path to the folder of interest
     * @return
     */
    public static List<String> listFolderContents(String folderName) {
        File[] contents = new File(folderName).listFiles();
        List<String> res = new ArrayList<>();
        for (File f : contents) {
            res.add(f.getAbsolutePath());
        }
        return (res);
    }

    /**
     * Loads a file and returns its contents as string
     *
     * @param filePath
     * @return
     * @throws IOException
     */
    public static String fileAsString(String filePath) throws IOException {
        StringBuilder result = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(filePath)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                result.append(line).append("\n");
            }
            return result.toString();
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves a string into a file
     *
     * @param toSave
     * @param filePath
     * @throws IOException
     */
    public static void saveStringToFile(String toSave, String filePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(toSave);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Performs file copy
     *
     * @param source
     * @param target
     * @throws IOException
     */
    public static void copyFile(String source, String target) throws IOException {
        Path sourcePath = Paths.get(source);
        Path targetPath = Paths.get(target);
        Files.copy(sourcePath, targetPath, REPLACE_EXISTING);
    }

    /**
     * For every label in a list, counts its occurrences in collection of
     * UniqueSequence objects. Full (nonunique) counts are returned
     *
     * @param sequences
     * @param labels
     * @return
     */
    private static Map<String, Integer> getTotalLabelCounts(Collection<UniqueSequence> sequences, List<String> labels) {
        Map<String, Integer> result = new HashMap<>();
        for (String label : labels) {
            result.put(label, 0);
        }
        for (UniqueSequence sequence : sequences) {
            for (Map.Entry<String, Integer> entry : sequence.getLabelsMap().entrySet()) {
                if (result.containsKey(entry.getKey())) {
                    int count = result.get(entry.getKey());
                    count += entry.getValue();
                    result.put(entry.getKey(), count);
                }
            }
        }
        return result;
    }

    /**
     * For every label in list, counts its occurrences in collection of
     * UniqueSequence objects. Each UniqueSequence having any positive count
     * associated with this label is counted as 1.
     *
     * @param sequences
     * @param labels
     * @return
     */
    private static Map<String, Integer> getUniqueLabelCounts(Collection<UniqueSequence> sequences, List<String> labels) {
        Map<String, Integer> result = new HashMap<>();
        for (String label : labels) {
            result.put(label, 0);
        }
        for (UniqueSequence sequence : sequences) {
            for (Map.Entry<String, Integer> entry : sequence.getLabelsMap().entrySet()) {
                if (result.containsKey(entry.getKey())) {
                    int count = result.get(entry.getKey());
                    count++;
                    result.put(entry.getKey(), count);
                }
            }
        }
        return result;
    }

    /**
     * Checks whether msa in "path" has enough conserved states satisfying
     * maximal gap proportion and minimal information content thresholds.
     *
     * @param alignmentLines clusters MSA
     * @param minMatchStates
     * @param minIc
     * @param maxGapProportion
     * @return
     * @throws IOException
     */
    public static boolean checkConservedStates(List<String> alignmentLines, int minMatchStates, double minIc, double maxGapProportion) throws IOException {
        List<Boolean> matchStates = defineMatchStates(alignmentLines, maxGapProportion, minIc, true); //true ensures conserved positions, not all match states between them are counted
        int count = 0;
        for (boolean position : matchStates) {
            if (position) {
                count++;
            }
        }
        return count >= minMatchStates;
    }

    /**
     * Returns list of information content values - one value for each column of
     * input alignment. Information contents are calculated based on
     * equiprobable model - no amino acid composition adjustments are performed.
     * Positions not having enough non-gap characters get information content of
     * -1.
     *
     * @param alignmentLines clusters MSA
     * @param maxGapProportion
     * @return
     * @throws IOException
     */
    public static List<Double> getInformationContents(List<String> alignmentLines, Double maxGapProportion) throws IOException {
        List<Map<Character, Integer>> positionLetterCounts = getPositionLetterCounts(alignmentLines);
        List<Double> result = new ArrayList<>();
        int seqCount = 0;
        for (Integer count : positionLetterCounts.get(0).values()) {
            seqCount += count;
        }

        for (Map<Character, Integer> positionCount : positionLetterCounts) {
            if ((positionCount.get('-') == null) || (((positionCount.get('-') + 0.0) / seqCount) <= maxGapProportion)) {
                result.add(getInformationContent(positionCount));
            } else {
                result.add(-1.0); //-1 for positions not having enough non-gap positions
            }
        }
        return result;
    }

    /**
     * Counts the occurrences of amino acids in each column of an MSA.
     *
     * @param alignmentLines A string containing a MSA.
     * @return A list, each member is a map respective to the appropriate
     * position in the MSA
     * @throws IOException
     */
    public static List<Map<Character, Integer>> getPositionLetterCounts(List<String> alignmentLines) throws IOException {
        List<Map<Character, Integer>> positionLetterCounts = new ArrayList<>();
        int lineLength = alignmentLines.get(0).trim().length();
        for (char aa : alignmentLines.get(0).trim().toCharArray()) {
            Map<Character, Integer> countMap = new HashMap<>();
            countMap.put(aa, 1);
            positionLetterCounts.add(countMap);
        }
        for (String line : alignmentLines.subList(1, alignmentLines.size())) {
            char[] aas = line.trim().toCharArray();
            for (int i = 0; i < lineLength; i++) {
                Integer count = positionLetterCounts.get(i).get(aas[i]);
                if (count == null) {
                    count = 0;
                }
                count++;
                positionLetterCounts.get(i).put(aas[i], count);
            }
        }
        return (positionLetterCounts);
    }

    /**
     * On the basis of criteria provided, decides which columns of a MSA fulfill
     * the criteria - maximal gap proportion and minimal information content and
     * will be considered match states.
     *
     * @param alignmentLines The MSA to be evaluated
     * @param maxGapProportion Maximal proportion of gap ('-' characters) in a
     * column to be considered a match state.
     * @param minIc Minimal information content of a MSA column to be considered
     * a match state
     * @param allowInnerGaps If true, a non-match state can lay between two
     * match states Otherwise, all the positions between the leftmost and the
     * rightmost match states are considered match states.
     * @return
     * @throws IOException
     */
    public static List<Boolean> defineMatchStates(List<String> alignmentLines, double maxGapProportion, double minIc, boolean allowInnerGaps) throws IOException {
        List<Double> informationContents = getInformationContents(alignmentLines, maxGapProportion);
        List<Boolean> result = new ArrayList<>();

        if (allowInnerGaps) {
            for (Double ic : informationContents) {
                if (ic >= minIc) {
                    result.add(true);
                } else {
                    result.add(false);
                }
            }
        } else {
            int minPosition = informationContents.size() + 1;
            int maxPosition = -1;
            for (int i = 0; i < informationContents.size(); i++) {
                if (informationContents.get(i) >= minIc) {
                    if (i < minPosition) {
                        minPosition = i;
                    }
                    if (i > maxPosition) {
                        maxPosition = i;
                    }
                }
            }
            for (int i = 0; i < informationContents.size(); i++) {
                if ((i >= minPosition) && (i <= maxPosition)) {
                    result.add(Boolean.TRUE);
                } else {
                    result.add(Boolean.FALSE);
                }
            }
        }

        return result;
    }

    /**
     * Transforms MSA in aln format to a2m format. Match columns (upper case)
     * are selected based on Global properties in Hammock. Saves the a2m file to
     * the temporal folder.
     *
     * @param cluster
     * @throws IOException
     * @throws InterruptedException
     */
    public static void alnToA2M(Cluster cluster) throws IOException, InterruptedException, DataException {
        aln2a2m(Settings.getInstance().getMsaDirectory() + cluster.getId() + ".aln", Settings.getInstance().getMsaDirectory() + cluster.getId() + ".a2m", Hammock.maxGapProportion, Hammock.minIc, Hammock.innerGapsAllowed);
    }

    private static void aln2a2m(String inFile, String outFile, Double maxGapProportion, double minIc, boolean allowInnerGaps) throws IOException, DataException {
        List<Boolean> matchStates = defineMatchStates(getAlignmentLines(inFile), maxGapProportion, minIc, allowInnerGaps);
        try (BufferedReader reader = new BufferedReader(new FileReader(inFile)); BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    writer.write(line + "\n");
                    line = reader.readLine();
                }
                if (line.length() != matchStates.size()) {
                    throw new DataException("Error. Wrong length of match state vector.");
                }
                StringBuilder newLine = new StringBuilder();
                for (int i = 0; i < line.length(); i++) {
                    if (line.charAt(i) == '-') {
                        if (matchStates.get(i)) {
                            newLine.append("-");
                        } else {
                            newLine.append(".");
                        }
                    } else {
                        if (matchStates.get(i)) {
                            newLine.append(Character.toUpperCase(line.charAt(i)));
                        } else {
                            newLine.append(Character.toLowerCase(line.charAt(i)));
                        }
                    }
                }
                writer.write(newLine.toString() + "\n");
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * transforms msa in selex format from inFile to msa in aln format and saves
     * it to outFile.
     *
     * @param inFile
     * @param outFile
     * @throws IOException
     */
    public static void selex2aln(String inFile, String outFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(inFile)); BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            reader.readLine(); //omit first line
            String line;
            while ((line = reader.readLine()) != null) {
                String[] splitLine = line.split("\\s+");
                StringBuilder newLine = new StringBuilder();
                for (char aa : splitLine[1].toCharArray()) {
                    switch (aa) {
                        case '.':
                            newLine.append('-');
                            break;
                        case '-':
                            newLine.append('-');
                            break;
                        default:
                            newLine.append(Character.toUpperCase(aa));
                            break;
                    }
                }
                writer.write(">" + splitLine[0] + "\n");
                writer.write(newLine.toString() + "\n");
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Gets the length of the first line of a file
     *
     * @param inFile
     * @return
     * @throws IOException
     */
    public static int getAlnLength(String inFile) throws IOException {
        String line;
        try (BufferedReader reader = new BufferedReader(new FileReader(inFile))) {
            reader.readLine();
            line = reader.readLine(); //second line
        } catch (IOException e) {
            throw new IOException(e);
        }
        return line.length();
    }

    /**
     * Checks whether alignment in inFile is longer or equal to maxLength
     *
     * @param alignmentLines clusters MSA
     * @param maxLength
     * @return
     * @throws IOException
     */
    public static boolean checkAlnLength(List<String> alignmentLines, int maxLength) throws IOException {
        return (alignmentLines.get(0).length() <= maxLength);
    }

    /**
     * Checks the number of inner gaps in the first and last lines of a msa
     *
     * @param alignmentLines
     * @param maxGaps
     * @return
     */
    public static boolean checkBothInnerGaps(List<String> alignmentLines, int maxGaps) {
        return ((countInnerGaps(alignmentLines.get(0)) <= maxGaps) && (countInnerGaps(alignmentLines.get(alignmentLines.size() - 1)) <= maxGaps));
    }

    private static int countInnerGaps(String line) {
        List<Integer> gapBlocks = new ArrayList<>();
        int currentBlock = 0;
        for (char aa : line.toCharArray()) {
            if (aa == '-') {
                currentBlock++;
            }
            if (aa != '-' && currentBlock > 0) {
                gapBlocks.add(currentBlock);
                currentBlock = 0;
            } //trailing gap block is never added
        }
        if ((line.toCharArray()[0]) == '-') { //starts with gaps
            gapBlocks = gapBlocks.subList(1, gapBlocks.size());
        }
        int sum = 0;
        for (int block : gapBlocks) {
            sum += block;
        }
        return sum;
    }

    /**
     * Calculates information content of one column in MSA. No amino acid
     * composition adjustments are performed.
     *
     * @param countMap
     * @return
     */
    private static double getInformationContent(Map<Character, Integer> countMap) {
        int size = 0;
        for (Map.Entry<Character, Integer> entry : countMap.entrySet()) {
            if (entry.getKey() != '-') {
                size += entry.getValue();
            }
        }
        List<Double> probabs = new ArrayList<>();
        for (Map.Entry<Character, Integer> entry : countMap.entrySet()) {
            if (entry.getKey() != '-') {
                probabs.add((entry.getValue() + 0.0) / size);
            }
        }
        double enthropy = 0.0;
        for (double probab : probabs) {
            enthropy += probab * (Math.log(probab) / Math.log(2)); //log base 2
        }
        return (-(Math.log(0.05) / Math.log(2)) + enthropy);
    }

    /**
     * Returns list of all UniqueSequence objects contained in a collection of
     * clusters
     *
     * @param clusters
     * @return
     */
    public static List<UniqueSequence> getAllSequences(Collection<Cluster> clusters) {
        List<UniqueSequence> result = new ArrayList<>();
        for (Cluster cl : clusters) {
            result.addAll(cl.getSequences());
        }
        return result;
    }
}

/**
 * Makes it possible to sort map according to values. Returns TOTAL ordering, so
 * it DOES NOT respect equals(). (it is not possible to respect equals in this
 * case - it would lead to key deletions)
 *
 * @author Adam Krejci
 */
class ValueComparator implements Comparator<String> {

    Map<String, Integer> base;

    public ValueComparator(Map<String, Integer> base) {
        this.base = base;
    }

    @Override
    public int compare(String a, String b) {
        if (base.get(a) >= base.get(b)) {
            return -1;
        } else {
            return 1;
        }
    }
}
