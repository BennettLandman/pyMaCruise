/**
 * Formula excerpt from VandyCRUISE MultipleToadSegmentation.java:
 * computeExactMemberships class-pair centroid/lesion similarity prelude.
 * This intentionally captures only the pairwise matrix construction.
 */
public final class ToadsPairSimilarityReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final byte OBJECT = 1, LESION = 4;
    private static final float PRIOR_SCALE = 0.17f;
    private static final float[][] CENTER = {
        {0.17f, 0.73f, 0.41f, 0.17f},
        {0.81f, 0.29f, 0.66f, 0.12f}
    };
    private static final float[][] STD = {
        {0.23f, 0.31f, 0.47f, 0.19f},
        {0.37f, 0.28f, 0.21f, 0.43f}
    };
    private static final float[] LESION_CENTER = {0.84f, 0.08f};
    private static final float[] LESION_STD = {0.13f, 0.34f};
    private static final float[] VARIANCE = {0.63f, 1.27f};
    private static final byte[] TYPES = {OBJECT, LESION, LESION, OBJECT};

    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static float[][] deviations(boolean zeroDeviation) {
        float[][] copy = new float[STD.length][];
        for (int c = 0; c < STD.length; c++) copy[c] = STD[c].clone();
        if (zeroDeviation) {
            copy[0][0] = 0.0f;
            copy[0][3] = 0.0f;
            copy[1][3] = 0.0f;
        }
        return copy;
    }

    private static void array(int[] values) {
        System.out.print('[');
        for (int i = 0; i < values.length; i++) {
            if (i != 0) System.out.print(',');
            System.out.print(values[i]);
        }
        System.out.print(']');
    }

    private static void floats(float[] values) {
        int[] output = new int[values.length];
        for (int i = 0; i < values.length; i++) output[i] = bits(values[i]);
        array(output);
    }

    private static String runCase(boolean zeroDeviation) {
        float[][] stddev = deviations(zeroDeviation);
        int n = TYPES.length, nc = CENTER.length;
        float[] similarity = new float[n*n];
        float[] centroidDistance = new float[n*n];
        float[] lesionDistance = new float[n*n];
        for (int k = 0; k < n; k++) for (int m = 0; m < n; m++) {
            int index = k*n+m;
            if (m != k) {
                float centroidDist, lesionDist = Float.NaN;
                if ((TYPES[k] == LESION) && (TYPES[m] == LESION)) centroidDist = 0.0f;
                else {
                    centroidDist = 0;
                    for (int c=0;c<nc;c++) centroidDist += VARIANCE[c]*(CENTER[c][k]-CENTER[c][m])*(CENTER[c][k]-CENTER[c][m])
                        /(stddev[c][k]+stddev[c][m])/(stddev[c][k]+stddev[c][m]);

                    if (TYPES[k] == LESION) {
                        lesionDist = 0;
                        for (int c=0;c<nc;c++) lesionDist += VARIANCE[c]*(LESION_CENTER[c]-CENTER[c][m])*(LESION_CENTER[c]-CENTER[c][m])
                            /(LESION_STD[c]+stddev[c][m])/(LESION_STD[c]+stddev[c][m]);
                    } else if (TYPES[m] == LESION) {
                        lesionDist = 0;
                        for (int c=0;c<nc;c++) lesionDist += VARIANCE[c]*(LESION_CENTER[c]-CENTER[c][k])*(LESION_CENTER[c]-CENTER[c][k])
                            /(LESION_STD[c]+stddev[c][k])/(LESION_STD[c]+stddev[c][k]);
                    }
                }
                similarity[index] = 1.0f/(1.0f + (1-PRIOR_SCALE)/PRIOR_SCALE*centroidDist);
                centroidDistance[index] = centroidDist;
                lesionDistance[index] = lesionDist;
            } else {
                similarity[index] = 0.0f;
                centroidDistance[index] = 0.0f;
                lesionDistance[index] = Float.NaN;
            }
        }
        StringBuilder out = new StringBuilder();
        out.append("{\"name\":\"").append(zeroDeviation ? "zero_deviation" : "asymmetric_mixed_pairs").append("\",");
        out.append("\"types\":[1,4,4,1],\"similarity_bits\":");
        int[] simBits = new int[similarity.length];
        int[] centroidBits = new int[centroidDistance.length];
        int[] lesionBits = new int[lesionDistance.length];
        for (int i=0;i<similarity.length;i++) {
            simBits[i]=bits(similarity[i]); centroidBits[i]=bits(centroidDistance[i]); lesionBits[i]=bits(lesionDistance[i]);
        }
        out.append(toJson(simBits)).append(",\"centroid_distance_bits\":").append(toJson(centroidBits));
        out.append(",\"lesion_distance_bits\":").append(toJson(lesionBits)).append('}');
        return out.toString();
    }

    private static String toJson(int[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i=0;i<values.length;i++) { if(i>0)out.append(',');out.append(values[i]); }
        return out.append(']').toString();
    }

    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships class-pair similarity prelude excerpt\",");
        System.out.print("\"source_sha256\":\""+SOURCE_SHA+"\",\"prior_scale_bits\":"+bits(PRIOR_SCALE)+",\"shape\":[4,4],\"cases\":[");
        System.out.print(runCase(false));
        System.out.print(',');
        System.out.print(runCase(true));
        System.out.println("]}");
    }
}
