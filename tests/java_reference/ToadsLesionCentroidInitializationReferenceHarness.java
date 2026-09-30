/*
 * Selected VandyCRUISE MultipleToadPreprocess.initialLesionCentroid excerpt,
 * as called by AlgorithmMultipleToads after segmentation initialization.
 * This validates scalar prior initialization only, not lesion segmentation.
 */
public final class ToadsLesionCentroidInitializationReferenceHarness {
    static final float[][] LOW = {{1.25f, 2.75f, 0.125f}, {0.001f, -1.25f, 20.0f}};
    static final float[][] HIGH = {{5.75f, 9.5f, 2.125f}, {0.005f, 1.25f, 20.75f}};
    static final float[][] PRIOR = {{0.2f, 0.7f, 0.9f}, {0.0f, 0.5f, 1.2f}};

    static float[] initialLesionCentroid(float[] low, float[] high, float[] prior) {
        float[] lesion = new float[low.length];
        for (int c=0; c<low.length; c++) {
            lesion[c] = low[c] / (high[c] - low[c]) + prior[c];
        }
        return lesion;
    }

    static int[] bits(float[] values) {
        int[] result = new int[values.length];
        for (int i=0; i<values.length; i++) result[i] = Float.floatToIntBits(values[i]);
        return result;
    }

    public static void main(String[] args) {
        System.out.println("{\"source_sha256\":{\"AlgorithmMultipleToads.java\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"MultipleToadPreprocess.java\":\"31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7\"},\"cases\":[");
        for (int i=0; i<LOW.length; i++) {
            float[] result = initialLesionCentroid(LOW[i], HIGH[i], PRIOR[i]);
            if (i>0) System.out.print(",");
            System.out.print("{\"low\":"+java.util.Arrays.toString(LOW[i])+",\"high\":"+java.util.Arrays.toString(HIGH[i])+",\"prior\":"+java.util.Arrays.toString(PRIOR[i])+",\"centroid_bits\":"+java.util.Arrays.toString(bits(result))+"}");
        }
        System.out.println("]}");
    }
}
