package hammock;

/**
 *
 * Class performs basic statistics computations
 */
public class Statistics {

    public static double mean(double[] vector) {
        double mean = 0;
        for (double num : vector) {
            mean += new Double(num);
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

    public static boolean checkCorrelation(Cluster cl1, Cluster cl2, double minCorrelation) {
        if ((Hammock.getLabels().size() < 2) || (minCorrelation == -1.0)) { //one label or any correlation => no control
            return true;
        }
        return (checkCorrelation(cl1.getLabelCountVector(), cl2.getLabelCountVector(), minCorrelation));
    }

    public static boolean checkCorrelation(Cluster cl, UniqueSequence seq, double minCorrelation) {
        if ((Hammock.getLabels().size() < 2) || (minCorrelation == -1)) { //one label or any correlation => no control
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

}


