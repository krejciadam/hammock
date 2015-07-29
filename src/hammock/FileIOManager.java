/*
 * Class handles most of work with files.
 */
package hammock;

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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Adam Krejci
 */
public class FileIOManager {

    /**
     * Returns a scoring matrix. Matrix uses standard amino acid ordering
     * (defined in UniqueSequence). Input file must contain scoring matrix in
     * BioJava format.
     *
     * @param matrixFilePath path to file with scoring matrix in BioJava format
     * @return two dimensional array representing scoring matrix.
     * @throws IOException
     * @throws HammockException
     */
    public static int[][] loadScoringMatrix(String matrixFilePath) throws IOException, HammockException {
        int[][] scoringMatrix = new int[24][24];

        Hammock.logger.logAndStderr("Loading scoring matrix...");

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
     * Loads unique sequences saved in a file in project's own format
     *
     * @param fileName path to file with unique sequences saved
     * @return List of UniqueSequence objects that were loaded from file.
     * @throws IOException
     * @throws FileFormatException
     */
    public static List<UniqueSequence> loadUniqueSequencesFromFile(String fileName) throws IOException, FileFormatException {
        List<UniqueSequence> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                result.add(processUniqueSequenceFileLine(line));
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new FileFormatException("Error in sequence file: " + fileName + ". File should contain"
                    + "lines in format: sequence|label1:count1|label2|count2|label3:count3| etc.");
        }
        return result;
    }

    /**
     * Loads a fasta file. Each seuqence must have a header line starting with
     * ">". Header may contain sequence label in format:
     * ">any_text|count|label". To each sequence not having label in the header,
     * default label "no_label" is assigned. Sequences in fasta file need not to
     * be unique.
     *
     * @param fileName Path to fasta file to be loaded
     * @return List of resulting unique sequences
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static List<UniqueSequence> loadUniqueSequencesFromFasta(String fileName) throws FileNotFoundException, IOException, FileFormatException {
        List<SeqWithLabel> tempList = new ArrayList<>();
        List<UniqueSequence> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line;
            String seq = null;
            String label = "no_label";
            int count = 1;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (seq != null) {
                        tempList.add(new SeqWithLabel(seq, label, count));
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
                    seq = line.trim();
                }
            }
            if (seq != null) {
                tempList.add(new SeqWithLabel(seq, label, count)); //add the last sequence
            }
        }
        Collections.sort(tempList); //alphabetic sort
        int startIndex = 0;
        Map<String, Integer> labelsMap = new HashMap<>();
        for (int i = 0; i < tempList.size(); i++) {
            SeqWithLabel currentSeq = tempList.get(i);
            if (!(currentSeq.getSeq().equals(tempList.get(startIndex).getSeq()))) {
                result.add(new UniqueSequence(tempList.get(startIndex).getSeq(), labelsMap));
                labelsMap = new HashMap<>();
                startIndex = i;
            }
            String lab = currentSeq.getLabel();
            Integer count = labelsMap.get(lab);
            if (count == null) {
                count = 0;
            }
            count += currentSeq.getCount();
            labelsMap.put(lab, count);
        }
        result.add(new UniqueSequence(tempList.get(startIndex).getSeq(), labelsMap)); //add last sequence
        return result;
    }

    /**
     * Loads UniqueSequence from a .csv table. Table must contain header
     * specifiing full list of labels. Table fields are separated by semicolon
     * (;)
     *
     * @param fileName
     * @param ignore
     * @return
     * @throws IOException
     */
    public static List<UniqueSequence> loadUniqueSequencesFromTable(String fileName, int[] ignore) throws IOException {
        Set<Integer> ignoreSet = new HashSet<>();
        List<UniqueSequence> result = new ArrayList<>();
        if (ignore != null) {
            for (int i : ignore) {
                ignoreSet.add(i);
            }
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line = reader.readLine(); //header line
            String[] headerLine = line.split(",");
            List<String> labels = new ArrayList<>();
            for (int i = 1; i < headerLine.length; i++) {
                if ((headerLine[i] != null) && (!(ignoreSet.contains(i)))) {
                    labels.add(headerLine[i]);
                }
            }
            while ((line = reader.readLine()) != null) {
                String[] splitLine = line.split(",");
                Map<String, Integer> labelsMap = new HashMap<>();
                int currentLabelsIndex = 0;
                for (int i = 1; i < splitLine.length; i++) {
                    if (!(ignoreSet.contains(i))) {
                        int value = Integer.decode(splitLine[i]);
                        if (value != 0) {
                            labelsMap.put(labels.get(currentLabelsIndex), value);
                        }
                        currentLabelsIndex++;
                    }
                }
                result.add(new UniqueSequence(splitLine[0], labelsMap));
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        return result;
    }

    /**
     * Legacy method, to be removed in next version
     *
     * @param fileName
     * @return
     * @throws IOException
     * @throws FileFormatException
     */
    public static List<Cluster> loadClustersFromFileOld(String fileName) throws IOException, FileFormatException {
        List<Cluster> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line;
            List<UniqueSequence> clusterSequences = null;
            Integer id = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if ((clusterSequences != null) && (id != null)) {
                        result.add(new Cluster(clusterSequences, id));
                    }
                    id = Integer.decode(line.substring(1));
                    clusterSequences = new ArrayList<>();
                } else {
                    clusterSequences.add(processUniqueSequenceFileLine(line));
                }
            }
            result.add(new Cluster(clusterSequences, id)); //add the last cluster
        } catch (NumberFormatException | NullPointerException e) {
            throw new FileFormatException("Error in cluster file: "
                    + fileName + " - wrong format. Original message: " + e.getMessage());
        }
        return result;
    }

    /**
     * Loads clusters from file. For every cluster, if MSA is saved, loads this
     * MSA as well and sets such cluster as hasMSA()
     *
     * @param fileName Path to file to be loaded
     * @return loaded clusters
     * @throws IOException
     * @throws FileFormatException
     */
    public static List<Cluster> loadClustersFromFile(String fileName) throws IOException, FileFormatException {
        List<Cluster> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String line;
            List<UniqueSequence> clusterSequences = null;
            StringBuilder alignment = null;
            Integer id = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("$")) {                  //read all alignment lines
                    alignment = new StringBuilder();
                    while (((line = reader.readLine()) != null) && (!(line.startsWith("&")))) {
                        alignment.append(line).append("\n");
                    }
                    alignment.deleteCharAt(alignment.length() - 1); //delete last /n
                    if (line == null) {
                        break;
                    }
                }
                if (line.startsWith("&")) {
                    if ((clusterSequences != null) && (id != null)) {
                        Cluster cl = new Cluster(clusterSequences, id);
                        if (alignment != null) {
                            FileIOManager.saveStringToFile(alignment.toString(), Settings.getInstance().getMsaDirectory() + id + ".aln");
                            alignment = null;
                            cl.setAsHasMSA();
                        }
                        result.add(cl);
                    }
                    id = Integer.decode(line.substring(1));
                    clusterSequences = new ArrayList<>();
                } else {
                    clusterSequences.add(processUniqueSequenceFileLine(line));
                }
            }
            Cluster cl = new Cluster(clusterSequences, id); //add the last cluster
            if (alignment != null) {
                FileIOManager.saveStringToFile(alignment.toString(), Settings.getInstance().getMsaDirectory() + id + ".aln");
                cl.setAsHasMSA();
            }
            result.add(cl);
        } catch (NumberFormatException | NullPointerException e) {
            e.printStackTrace();
            throw new FileFormatException("Error in cluster file: "
                    + fileName + " - wrong format. Original message: " + e.getMessage());
        }
        return result;
    }

    /**
     * Processes one line of a file with UniqueSequence objects saved
     *
     * @param line line from file with UniqueSequence objects saved
     * @return
     * @throws ArrayIndexOutOfBoundsException
     */
    private static UniqueSequence processUniqueSequenceFileLine(String line) throws ArrayIndexOutOfBoundsException {
        Map<String, Integer> labelMap = new HashMap<>();
        String[] splitLine = line.split("\\|");
        String sequence = splitLine[0];
        for (int i = 1; i < splitLine.length; i++) {
            String[] splitField = splitLine[i].split(":");
            labelMap.put(splitField[0], Integer.parseInt(splitField[1]));
        }
        return new UniqueSequence(sequence, labelMap);
    }

    /**
     * Saves a csv file containing one line per sequence for every cluster on
     * input.
     *
     * @param clusters Clusters to be saved
     * @param filePath Path to resulting file
     * @param labels List of labels. Only these labels and their appropriate
     * counts will be saved in resulting file.
     * @throws IOException
     */
    public static void saveClusterSequencesToCsv(Collection<Cluster> clusters, String filePath, List<String> labels) throws IOException {
        List<Cluster> sortedList = new ArrayList<>(clusters);
        Collections.sort(sortedList, Collections.reverseOrder());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("cluster_id;sequence;alignment;sum");
            for (String label : labels) {
                writer.write(";" + label);
            }
            writer.newLine();
            for (Cluster cl : sortedList) {
                List<String> msa = null;
                Map<String, String> msaMap = null;
                if (cl.hasMSA()) {
                    msa = FileIOManager.getAlignmentLines(cl);
                    msaMap = new HashMap<>();
                    for (String msaLine : msa) {
                        msaMap.put(msaLine.replace("-", ""), msaLine);
                    }
                }
                List<UniqueSequence> sequences = cl.getSequences();
                Collections.sort(sequences, Collections.reverseOrder(new UniqueSequenceSizeAlphabeticComparator()));
                for (int i = 0; i < sequences.size(); i++) {
                    UniqueSequence seq = sequences.get(i);
                    writer.write(cl.getId() + ";" + seq.getSequenceString() + ";");
                    if (cl.hasMSA()) {
                        writer.write(msaMap.get(seq.getSequenceString()) + ";");
                    } else {
                        writer.write("NA;");
                    }
                    writer.write("" + seq.size());
                    for (String label : labels) {
                        Integer count = seq.getLabelsMap().get(label);
                        if (count == null) {
                            count = 0;
                        }
                        writer.write(";" + count);
                    }
                    writer.newLine();
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves cluster overview to .csv file
     *
     * @param clusters
     * @param filePath
     * @param labels
     * @throws IOException
     */
    public static void SaveClustersToCsv(Collection<Cluster> clusters, String filePath, List<String> labels) throws IOException {
        List<Cluster> sortedList = new ArrayList<>(clusters);
        Collections.sort(sortedList, Collections.reverseOrder());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("cluster_id;main_sequence;sum");
            for (String label : labels) {
                writer.write(";" + label);
            }
            writer.newLine();
            for (Cluster cl : sortedList) {
                List<UniqueSequence> sequences = cl.getSequences();
                Collections.sort(sequences, Collections.reverseOrder());
                writer.write(cl.getId() + ";" + sequences.get(0).getSequenceString() + ";" + cl.size());
                Map<String, Integer> clusterCountMap = getClusterLabelsMap(cl);
                for (String label : labels) {
                    Integer count = clusterCountMap.get(label);
                    if (count == null) {
                        count = 0;
                    }
                    writer.write(";" + count);
                }
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException(e);
        }

    }

    /**
     * Returns labels count map for whole cluster (sum of all sequences)
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
     * Saves file containing input statistics - information on sequences and
     * their counts
     *
     * @param sequences
     * @param labels
     * @param filePath
     * @throws IOException
     */
    public static void saveInputStatistics(Collection<UniqueSequence> sequences, List<String> labels, String filePath) throws IOException {
        Map<String, Integer> totalMap = getTotalLabelCounts(sequences, labels);
        Map<String, Integer> uniqueMap = getUniqueLabelCounts(sequences, labels);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (String label : labels) {
                writer.write(";" + label);
            }
            writer.newLine();
            writer.write("total_count");
            for (String label : labels) {
                writer.write(";" + totalMap.get(label));
            }
            writer.newLine();
            writer.write("unique_count");
            for (String label : labels) {
                writer.write(";" + uniqueMap.get(label));
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    /**
     * Saves a Collection of Savable objects to a file in reversed natural
     * ordering. No linebreaks are added.
     *
     * @param clusters
     * @param filePath Path to a file Savable objects will be saved to.
     * @throws IOException
     */
    public static void saveClustersToFile(Collection<Cluster> clusters, String filePath) throws IOException {
        List<Cluster> sortedList = new ArrayList<>(clusters);
        Collections.sort(sortedList, Collections.reverseOrder());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (Cluster cluster : sortedList) {
                writer.write(cluster.getSavableString());
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

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
     * Loads temporal MSA file for a cluster and returns onlu lines containing
     * alignment in form of a (multiline) string
     *
     * @param cl
     * @return
     * @throws IOException
     */
    private static List<String> getAlignmentLines(Cluster cl) throws IOException {
        return FileIOManager.getAlignmentLines(Settings.getInstance().getMsaDirectory() + cl.getId() + ".aln");
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
     * Saves one Cluster object to a fasta file. No linebreaks are addeded.
     * Cluster's sequences are saved in reversed natural order.
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
     * Merges MSA files for two clusters, adds any gaps specified
     *
     * @param cl1 first cluster to be merged
     * @param cl2 second cluster to be merged
     * @param gaps1 gaps to be added into first cluster
     * @param gaps2 gaps to be added into second cluster
     * @param newId id of resulting cluster
     * @throws IOException
     */
    public static void mergeAlignedClusters(Cluster cl1, Cluster cl2, List<Integer> gaps1, List<Integer> gaps2, int newId) throws IOException {
        String cl1String = insertGapsIntoAlignment(cl1.getId(), gaps1, newId, 1);
        String cl2String = insertGapsIntoAlignment(cl2.getId(), gaps2, newId, cl1.getUniqueSize() + 1);
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
    private static String insertGapsIntoAlignment(int id, List<Integer> gapPositions, int newId, int startingSequenceId) throws IOException {
        StringBuilder res = new StringBuilder();
        int currentId = startingSequenceId;
        String[] alnLines = fileAsString(Settings.getInstance().getMsaDirectory() + id + ".aln").split("\n");
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

    public static void makeShiftedClusterAlignment(UniqueSequence pivot, List<AligningScorerResult> results, int clusterId) throws IOException {
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
        saveStringToFile(alignment.toString(), Settings.getInstance().getMsaDirectory() + clusterId + ".aln");
    }

    /**
     * For every cluster, saves path to its appropriate .hmm file in temporal
     * files
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
     * Writes a collection of sequences to a fasta file. Assigns a unique
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
     * >newId_n where n is sequence index (from top to botton in file)
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
     * Deletes all files in specified folder
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
     * @param folderName path to the folder of interest
     * @return 
     */
    public static List<String> listFolderContents(String folderName){
        File[] contents = new File(folderName).listFiles();
        List<String> res = new ArrayList<>();
        for (File f : contents){
            res.add(f.getAbsolutePath());
        }
        return(res);
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
     * For every label in list, counts its occurrences in collection of
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
     * Checks whether msa in "path" has enough match states satisfying maximal
     * gap proportion and minimal information content thresholds.
     *
     * @param alignmentLines clusters MSA
     * @param minMatchStates
     * @param minIc
     * @param maxGapProportion
     * @return
     * @throws IOException
     */
    public static boolean checkMatchStatesAndIc(List<String> alignmentLines, int minMatchStates, double minIc, double maxGapProportion) throws IOException {
        List<Boolean> matchStates = defineMatchStates(alignmentLines, maxGapProportion, minIc);
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

    public static List<Boolean> defineMatchStates(List<String> alignmentLines, double maxGapProportion, double minIc) throws IOException {
        List<Double> informationContents = getInformationContents(alignmentLines, maxGapProportion);
        List<Boolean> result = new ArrayList<>();
        for (Double ic : informationContents) {
            if (ic >= minIc) {
                result.add(true);
            } else {
                result.add(false);
            }
        }
        return result;
    }

    /**
     * Transforms MSA in aln format to a2m format. Match columns (upper case)
     * are selected based on maximal proportion of gaps and minimal information
     * content in a match state.
     *
     * @param inFile
     * @param outFile
     * @param maxGapProportion
     * @param minIc
     * @throws IOException
     * @throws DataException
     */
    public static void aln2a2m(String inFile, String outFile, Double maxGapProportion, double minIc) throws IOException, DataException {
        List<Boolean> matchStates = defineMatchStates(getAlignmentLines(inFile), maxGapProportion, minIc);
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

    public static boolean checkBothInnerGaps(List<String> alignmentLines, int maxGaps){
        return((countInnerGaps(alignmentLines.get(0)) <= maxGaps) && (countInnerGaps(alignmentLines.get(alignmentLines.size() - 1)) <= maxGaps));
    }

    public static boolean checkLastInnerGaps(String inFile, int maxGaps) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(inFile))) {
            String lastLine = null;
            String line;
            while ((line = reader.readLine()) != null) {
                lastLine = line;
            }
            return (countInnerGaps(lastLine) <= maxGaps);

        } catch (IOException e) {
            throw new IOException(e);
        }
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
 * This class is a pre-stage for UniqueSequence object when loading sequences
 * from .fasta files
 *
 * @author Adam Krejci
 */
class SeqWithLabel implements Comparable<SeqWithLabel> {

    private final String seq;
    private final String label;
    private final int count;

    public SeqWithLabel(String seq, String label, int count) {
        this.seq = seq;
        this.label = label;
        this.count = count;
    }

    public String getSeq() {
        return seq;
    }

    public String getLabel() {
        return label;
    }

    public int getCount() {
        return count;
    }

    @Override
    public int compareTo(SeqWithLabel o) {
        return this.seq.compareTo(o.getSeq());
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
