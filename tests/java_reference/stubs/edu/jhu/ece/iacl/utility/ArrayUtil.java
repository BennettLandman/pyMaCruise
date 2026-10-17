package edu.jhu.ece.iacl.utility;

/** Small headless subset used by the selected AtlasEM rigid source. */
public final class ArrayUtil {
    private ArrayUtil() {}

    public static float[] toFloat(double[] values) {
        float[] result = new float[values.length];
        for (int i = 0; i < values.length; i++) result[i] = (float) values[i];
        return result;
    }

    public static double[][] toDouble(float[][] values) {
        double[][] result = new double[values.length][];
        for (int i = 0; i < values.length; i++) {
            result[i] = new double[values[i].length];
            for (int j = 0; j < values[i].length; j++) result[i][j] = values[i][j];
        }
        return result;
    }
}
