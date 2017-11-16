package hammock;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * Class performs statistics computations
 */
public class Statistics {

    private static final double BETA = 200;
    private static final double SIGMA = 10;
    private static final double MATRIX_SCALE_FACTOR = 2.88539;

    private static final double[] backgroundProbabilities = new double[]{0.074, 0.052, 0.045, 0.054, 0.025,
        0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058,
        0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032,
        0.073};

    private static final Map<Character, Double> background;
    private static Map<Character, Map<Character, Double>> frequencyMatrix = null;

    static {
        Map<Character, Double> BACKGROUND = new HashMap<>();
        for (int i = 0; i < UniqueSequence.getAminoAcids().length; i++) {
            BACKGROUND.put(UniqueSequence.getAminoAcids()[i], backgroundProbabilities[i]);
        }
        background = Collections.unmodifiableMap(BACKGROUND);
        try {
            frequencyMatrix = FileIOManager.loadFrequencyMatrix(Hammock.countMatrixFile);
        } catch (IOException ex) {
            Hammock.logger.logAndStderr("Error. Can't load frequency matrix");
            System.exit(1);
        }
    }

    public static double mean(double[] vector) {
        double mean = 0;
        for (double num : vector) {
            mean += num;
        }
        mean /= vector.length;
        return mean;
    }

    public static double pearsonCorrelation(double[] vector1, double[] vector2) {
        double mean1 = mean(vector1);
        double mean2 = mean(vector2);
        double top = 0, down1 = 0, down2 = 0;
        for (int i = 0; i < vector1.length; i++) {
            double a = vector1[i] - mean1;
            double b = vector2[i] - mean2;
            top += (a * b);
            down1 += a * a;
            down2 += b * b;
        }
        double r = top / (Math.sqrt(down1) * Math.sqrt(down2));
        return r;
    }

    public static double pearsonCorrelation(int[] vector1, int[] vector2) {
        return (pearsonCorrelation(intToDoubleArray(vector1), intToDoubleArray(vector2)));
    }

    private static double[] intToDoubleArray(int[] source) {
        double[] res = new double[source.length];
        for (int i = 0; i < source.length; i++) {
            res[i] = source[i];
        }
        return res;
    }
    
    /**
     * Checks the the correlations of the label vectors of two clusters
     */
    public static boolean checkCorrelation(Cluster cl1, Cluster cl2, double minCorrelation) {
        if ((Hammock.getLabels().size() < 2) || (minCorrelation == -1.0)) { //one label or any correlation => no control
            return true;
        }
        return (checkCorrelation(cl1.getLabelCountVector(), cl2.getLabelCountVector(), minCorrelation));
    }

    /**
     * Checks the correlations of label vectors of a cluster and a sequence
     */
    public static boolean checkCorrelation(Cluster cl, UniqueSequence seq, double minCorrelation) {
        if ((Hammock.getLabels().size() < 2) || (minCorrelation <= -1.0)) { //one label or any correlation => no control
            return true;
        }
        return (checkCorrelation(cl.getLabelCountVector(), getLabelsVector(seq), minCorrelation));
    }

    private static boolean checkCorrelation(int[] vec1, int[] vec2, double minCorrelation) {
        if ((Hammock.getLabels().size() < 2) || (minCorrelation == -1)) { //one label or any correlation => no control
            return true;
        }
        double r = Statistics.pearsonCorrelation(vec1, vec2);
        if (r == Double.NaN) { //Clusters with undefined correlation can be merged
            return true;
        }
        return (r >= minCorrelation);
    }

    private static int[] getLabelsVector(UniqueSequence seq) {
        int[] seqLabelCounts = new int[Hammock.getLabels().size()];
        for (int i = 0; i < Hammock.getLabels().size(); i++) {
            String label = Hammock.getLabels().get(i);
            Integer value = seq.getLabelsMap().get(label);
            if (value == null) {
                value = 0;
            }
            seqLabelCounts[i] = value;
        }
        return seqLabelCounts;
    }
    
    /**
     * Not in use now 
     */
    public static boolean checkPerMatchStateKldIncrement(List<String> newClusterLines, List<String> cl1Lines, List<String> cl2Lines) throws IOException, DataException{
        List<Boolean> newClusterMatchStates = FileIOManager.defineMatchStates(newClusterLines, Hammock.maxGapProportion, Hammock.minIc, Hammock.innerGapsAllowed);
        List<Boolean> cl1MatchStates = FileIOManager.defineMatchStates(cl1Lines, Hammock.maxGapProportion, Hammock.minIc, Hammock.innerGapsAllowed);
        List<Boolean> cl2MatchStates = FileIOManager.defineMatchStates(cl2Lines, Hammock.maxGapProportion, Hammock.minIc, Hammock.innerGapsAllowed);
        List<Double> newClusterKlds = getClusterKlds(newClusterLines, newClusterMatchStates);
        List<Double> cl1Klds = getClusterKlds(cl1Lines, cl1MatchStates);
        List<Double> cl2Klds = getClusterKlds(cl2Lines, cl2MatchStates);
        int oldCount = cl1Klds.size() + cl2Klds.size();
        double cl1KldSum = sum(cl1Klds);
        int cl1MatchStatesCount = 0;
        for (Boolean b : cl1MatchStates){
            if (b){
                cl1MatchStatesCount++;
            }
        }
        cl1KldSum = cl1KldSum / (cl1MatchStatesCount + 0.0);
        double cl2KldSum = sum(cl2Klds);
        int cl2MatchStatesCount = 0;
        for (Boolean b : cl2MatchStates){
            if (b){
                cl2MatchStatesCount++;
            }
        }
        cl2KldSum = cl2KldSum / (cl2MatchStatesCount + 0.0);
        double oldAverage = (cl1KldSum + cl2KldSum) / (oldCount + 0.0);
        
        double newAverage = sum(newClusterKlds);
        int newCount = newClusterKlds.size();
        int newClusterMatchStatesCount = 0;
        for (Boolean b : newClusterMatchStates){
            if (b){
                newClusterMatchStatesCount++;
            }
        }
        newAverage = newAverage / (newClusterMatchStatesCount + 0.0);
        newAverage = newAverage / (newCount + 0.0);
        return(newAverage >= oldAverage);
    }
    
    
    /**
     * Returns the average KLD for a system of clusters
     * @param alnFiles paths to cluster MSAs
     * @param allPositions all positions are counted. If false, only match states are counted
     * @return
     * @throws IOException
     * @throws DataException 
     */
    public static double getMeanSystemKld(Collection<String> alnFiles, boolean allPositions) throws IOException, DataException{
        List<Double> klds = new ArrayList<>();
        for (String file : alnFiles){
            klds.addAll(getClusterKlds(file, allPositions));
        }
        return(sum(klds) / (klds.size() + 0.0));
    }
    
    /**
     * Returns the average KLD of a cluster
     * @param alignmentPath path to clusters MSA
     * @param allPositions all positions are counted. If false, only match states are counted
     * @return
     * @throws DataException
     * @throws IOException 
     */
    public static double getMeanClusterKld(String alignmentPath, boolean allPositions) throws DataException, IOException{
        List<Double> klds = getClusterKlds(alignmentPath, allPositions);
        return(sum(klds) / (klds.size() + 0.0));
    }
    
    /**
     * Returns a list of KLDs for a cluster, one for each peptide in the cluster.
     * @param alignmentPath path to clusters MSA
     * @param allPositions all positions are counted. If false, only match states are counted
     * @return
     * @throws IOException
     * @throws DataException 
     */
    public static List<Double> getClusterKlds(String alignmentPath, boolean allPositions) throws IOException, DataException{
        List<String> alignmentLines = FileIOManager.getAlignmentLines(alignmentPath);
        List<Boolean> matchStates = new ArrayList<>();
        if (allPositions){
            for (int i = 0; i < alignmentLines.get(0).length(); i++){
                matchStates.add(Boolean.TRUE);
            }
        } else{
            matchStates = FileIOManager.defineMatchStates(alignmentLines, Hammock.maxGapProportion, Hammock.minIc, Hammock.innerGapsAllowed);
        }
        return(getClusterKlds(alignmentLines, matchStates));
    }
    
    private static List<Double> getClusterKlds(List<String> alignmentLines, List<Boolean> matchStates) throws IOException, DataException{
        List<Double> result = new ArrayList<>();
        List<Map<Character, Integer>> letterCounts = FileIOManager.getPositionLetterCounts(alignmentLines);
        for (String line : alignmentLines){
            result.add(getPeptideKld(letterCounts, line, matchStates));
        }
        return(result);
    }

    /**
     * 
     * @param letterCounts Cluster position-wise letter counts INCLUDING the alignedPeptide
     * @param alignedPeptide The line to be scored. Must originate from the cluster used to calculate letterCounts
     * @param matchStates
     * @return
     * @throws IOException
     * @throws DataException 
     */
    private static Double getPeptideKld(List<Map<Character, Integer>> letterCounts, String alignedPeptide, List<Boolean> matchStates) throws IOException, DataException {
        if (alignedPeptide.length() != matchStates.size()) {
            throw new DataException("Error. Match states vector and aligned peptide lengths differ.");
        }
        
        Double result = 0.0;
        for (int i = 0; i < letterCounts.size(); i++) {
            Double positionKld;
            if (matchStates.get(i)) {
                Map<Character, Integer> countMap = new HashMap<>(letterCounts.get(i));
                Character peptideAA = alignedPeptide.charAt(i);
                if (peptideAA == '-') {
                    result += 0.0;
                    continue;
                }
                if (countMap.containsKey(peptideAA) && (countMap.get(peptideAA) != 0)) { //for intra-cluster comparisons, this should always be true (always at least 1)
                    countMap.put(peptideAA, (countMap.get(peptideAA) - 1));
                } else {
                    throw new DataException("Error. Scored peptide probably not from scored cluster");
                }
                double sequenceCount = 0;
                for (int c : countMap.values()) {
                    sequenceCount += c;
                }
                if (countMap.containsKey('-') && (sequenceCount == countMap.get('-'))){
                    result += 0.0;
                    continue;
                }
                Map<Character, Double> correctedFrequencies = getCorrectedFrequencies(countMap);
                positionKld = Math.log(correctedFrequencies.get(peptideAA) / background.get(peptideAA)) * ((sequenceCount) / (sequenceCount + SIGMA));
                positionKld = positionKld * MATRIX_SCALE_FACTOR;
                result += positionKld;
            }
        }
        return (result);
    }

    private static Map<Character, Double> getCorrectedFrequencies(Map<Character, Integer> counts) {
        Map<Character, Double> result = new HashMap<>();
        double sum = 0.0;
        for (int i : counts.values()) {
            sum += i;
        }
        for (Character aa : UniqueSequence.getAminoAcids()) {
            double gi = 0;
            for (Map.Entry<Character, Integer> entry : counts.entrySet()) {
                if (entry.getKey() != '-') { //gaps add no weight
                    double fj = (entry.getValue() + 0.0) / sum;
                    double qij = frequencyMatrix.get(entry.getKey()).get(aa);
                    gi += fj  * qij;
                }
            }
            double fi;
            if (counts.containsKey(aa)) {
                fi = (counts.get(aa) + 0.0) / sum;
            } else {
                fi = 0.0;
            }
            double Qi = (((sum - 1) * fi) + (BETA * gi)) / ((sum - 1) + BETA);
            result.put(aa, Qi);
        }
        return (result);
    }
    
    /** NOT in use now
     * For all the clusters in the folder of cluster MSAs, checks which sequence belongs to 
     * which cluster.
     */
    public static List<Integer> getOrderedClusterMemberships(String folder) throws IOException{
        List<String> sequences = new ArrayList<>();
        Map<String, Integer> membershipMap = new HashMap<>();
        int id = 0;
        for (String file : FileIOManager.listFolderContents(folder)){
            for (String sequence : FileIOManager.getAlignmentLines(file)){
                sequence = sequence.replace("-", "");
                sequences.add(sequence);
                membershipMap.put(sequence, id);
            }
            id++;
        }
        Collections.sort(sequences);
        List<Integer> result = new ArrayList<>();
        for(String s : sequences){
            result.add(membershipMap.get(s));
        }
        return(result);
    }
    
    /**
     * Returns sum of a list of doubles
     * @param list the list to be summed
     * @return 
     */
    public static double sum(Collection<Double> list){
        double result = 0.0;
        for (Double d : list){
            result += d;
        }
        return(result);
    }
    
    /**
     * Returns a list of values representing the x-axis.
     * Sorts from largest to smallest
     * @param scores
     * @return 
     */
    private static List<Double> getXAxis(List<Double> scores){
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (Double d : scores){
            stats.addValue(d);
        }
        List<Double> result = new ArrayList<>();
        long min = Math.round(stats.getMin() * 10.0);
        long max = Math.round(stats.getMax() * 10.0);
        for (long i = max; i >= min; i -= 1){
            result.add((i + 0.0) / 10);
        }
        return(result);
    }
    
    /**
     * Histogram
     * @param axis
     * @param scores
     * @return 
     */
    private static List<Double> toHist(List<Double> axis, List<Double> scores){
        List<Double> result = new ArrayList<>();
        for (Double x : axis){
            result.add(Collections.frequency(scores, x) + 0.0);
        }
        return(result);
    }
    
    /**
     * Simple moving-average smoothing
     * @param vector
     * @param bandwidth
     * @return 
     */
    private static List<Double> smooth(List<Double> vector, int bandwidth){
        List<Double> result = new ArrayList<>();
        for (int i = 0; i <= (vector.size() - bandwidth); i++){
            result.add(sum(vector.subList(i, i + bandwidth)) / bandwidth);
        }
        return(result);
    }
    
    /**
     * Returns the automatically derived threshold
     * @param scores
     * @param bandwidth
     * @param tolerance
     * @param leaveout
     * @return 
     */
    public static double getThreshold(List<Double> scores, int bandwidth, double tolerance, int leaveout){
        List<Double> axis = getXAxis(scores);
        List<Double> smoothed = smooth(toHist(axis, scores), bandwidth);
        List<Double> unitAxis = new ArrayList<>();
        for (int i = 0; i < smoothed.size(); i++){
            unitAxis.add(i + 0.0);
        }
        List<Double> unitAxisLeaveout = leaveout(unitAxis, leaveout);
        List<Double> smoothedLeaveout = leaveout(smoothed, leaveout);
        PolynomialSplineFunction interpolated = interpolate(unitAxisLeaveout, smoothedLeaveout);
        PolynomialSplineFunction derivative = interpolated.polynomialSplineDerivative();
        PolynomialSplineFunction derivative2 = derivative.polynomialSplineDerivative();
        int index = 0;
        for (double i = 0; i < unitAxisLeaveout.get(unitAxisLeaveout.size() - 1); i += 0.01){
            double val = derivative.value(i);
            if (Math.abs(val) <= 0.0025){ //local exterm
//                System.out.println(i + " " + val);
                double val2 = derivative2.value(i);
//                System.out.println(val2);
                if (val2 > 0.0){ //local minimum
                    boolean accept = false;
                    for (double j = 0.01; j <= 10 && j + i < unitAxisLeaveout.get(unitAxisLeaveout.size() - 1); j += 0.01){
                        if (derivative.value(i + j) >= tolerance){
                            accept = true;
                        }
                    }
                    if (accept){
                        index = (int)Math.round(i);
                        break;
                    }
                }
            }
        }
        return(axis.get(index));
    }
    
    /**
     * Takes only each ?-th member of the list
     * @param list
     * @param leaveout
     * @return 
     */
    private static List<Double> leaveout(List<Double> list, int leaveout){
        List<Double> res = new ArrayList<>();
        int left = Integer.MAX_VALUE;
        for (int i = 0; i < list.size(); i++){
            if (left > leaveout){
                res.add(list.get(i));
                left = 0;
            }
            left++;
        }
        return(res);
    }
                

    /**
     * Returns an interpolation spline.
     * @param axis
     * @param values
     * @return 
     */
    private static PolynomialSplineFunction interpolate(List<Double> axis, List<Double> values){
        SplineInterpolator interpolator = new SplineInterpolator();
        Double[] x = new Double[values.size()];
        x = axis.subList(0, values.size()).toArray(x);
        Double[] y = new Double[values.size()];
        y = values.toArray(y);
        PolynomialSplineFunction interpolated = interpolator.interpolate(ArrayUtils.toPrimitive(x), ArrayUtils.toPrimitive(y));
        return(interpolated);
    }    
}