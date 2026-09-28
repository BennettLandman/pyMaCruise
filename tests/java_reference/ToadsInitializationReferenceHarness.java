/*
 * Test-only copy of selected MultipleToadPreprocess initialization formulas.
 *
 * Source: toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadPreprocess.java,
 * methods initialCentroids and initialExactMemberships. The selected
 * AlgorithmMultipleToads constructor sets useApprox=false. This harness keeps
 * Java float arithmetic and iteration order and has no MIPAV/JIST dependency.
 */
import java.util.Locale;

public final class ToadsInitializationReferenceHarness {
    private static final int CHANNELS = 2;
    private static final int NX = 6, NY = 6, NZ = 6, CLASSES = 4;
    private static final int BORDER = 1;
    private static final float ZERO = 1e-20f;
    private static final float INF = 1e20f;
    private static final float[] LOW = {2.0f, 3.0f};
    private static final float[] HIGH = {10.0f, 13.0f};
    private static final float[][] PRIORS = {
        {0.1f, 0.4f, 0.4f, 0.8f},
        {0.2f, 0.5f, 0.7f, 0.9f}
    };

    private static float[][][][] images(float[][] centroids) {
        float[][][][] result = new float[CHANNELS][NX][NY][NZ];
        for (int c = 0; c < CHANNELS; c++) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int z = 0; z < NZ; z++) {
                        int code = c == 0
                            ? (x * 11 + y * 7 + z * 3) % 23
                            : (x * 5 + y * 13 + z * 17) % 29;
                        result[c][x][y][z] = 0.2f + code / (c == 0 ? 20.0f : 25.0f);
                    }
                }
            }
        }
        // A zero-distance sample exercises the source 1e20 reciprocal cap.
        result[0][2][2][2] = centroids[0][1];
        result[1][2][2][2] = centroids[1][1];
        return result;
    }

    private static float[][] initialCentroids() {
        float[][] result = new float[CHANNELS][CLASSES];
        for (int c = 0; c < CHANNELS; c++) {
            float span = HIGH[c] - LOW[c];
            float offset = LOW[c] / span;
            for (int k = 0; k < CLASSES; k++) result[c][k] = offset + PRIORS[c][k];
        }
        return result;
    }

    private static float[][][][] initialMemberships(float[][][][] image, float[][] centroid) {
        float[][][][] mems = new float[NX][NY][NZ][CLASSES];
        float[] grouping = new float[CLASSES];
        for (int k = 0; k < CLASSES; k++) {
            grouping[k] = 0.0f;
            for (int m = 0; m < CLASSES; m++) {
                if (centroid[0][m] == centroid[0][k]) grouping[k]++;
            }
        }
        for (int x = BORDER; x < NX - BORDER - 1; x++) {
            for (int y = BORDER; y < NY - BORDER - 1; y++) {
                for (int z = BORDER; z < NZ - BORDER - 1; z++) {
                    float den = 0.0f;
                    float[] num = new float[CLASSES];
                    for (int k = 0; k < CLASSES; k++) {
                        num[k] = 0.0f;
                        for (int c = 0; c < CHANNELS; c++) {
                            float diff = image[c][x][y][z] - centroid[c][k];
                            num[k] += diff * diff;
                        }
                        if (num[k] > ZERO) num[k] = 1.0f / num[k];
                        else num[k] = INF;
                        den += num[k];
                    }
                    if (den > 0.0f) {
                        for (int k = 0; k < CLASSES; k++) {
                            mems[x][y][z][k] = grouping[k] * num[k] / den;
                        }
                    }
                }
            }
        }
        return mems;
    }

    private static void append(StringBuilder output, float value) {
        output.append(Float.toString(value));
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float[][] centroid = initialCentroids();
        float[][][][] image = images(centroid);
        float[][][][] membership = initialMemberships(image, centroid);
        StringBuilder out = new StringBuilder();
        out.append("{\"shape\":[").append(NX).append(',').append(NY).append(',').append(NZ)
            .append("],\"crop_border\":").append(BORDER).append(",\"centroids\":[");
        for (int c = 0; c < CHANNELS; c++) {
            if (c > 0) out.append(',');
            out.append('[');
            for (int k = 0; k < CLASSES; k++) {
                if (k > 0) out.append(',');
                append(out, centroid[c][k]);
            }
            out.append(']');
        }
        out.append("],\"memberships\":[");
        boolean first = true;
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            for (int k = 0; k < CLASSES; k++) {
                if (!first) out.append(',');
                append(out, membership[x][y][z][k]);
                first = false;
            }
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
