/** Selected-method excerpt of VandyCRUISE MultipleToadSegmentation.generateClassification. */
public final class ToadsGenerateClassificationReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX = 3, NY = 2, NZ = 4, CLASSES = 4;
    private static final byte[] TEMPLATE = new byte[] {-128, -1, 42, 127};

    private static final float[][][][] MEMS = new float[NX][NY][NZ][CLASSES];
    private static final int[][][] BEST = new int[NX][NY][NZ];

    private static void set(int linearXFast, int k, float value) {
        int x = linearXFast % NX;
        int y = (linearXFast / NX) % NY;
        int z = linearXFast / (NX * NY);
        MEMS[x][y][z][k] = value;
    }

    private static void initialize() {
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            for (int k = 0; k < CLASSES; k++) MEMS[x][y][z][k] = -1.0f;
            BEST[x][y][z] = (x + 2 * y + 3 * z) % CLASSES;
        }
        // All-nonpositive, NaN-only, +infinity, strict tie, and class-order cases.
        for (int k = 0; k < CLASSES; k++) set(0, k, 0.0f);
        for (int k = 0; k < CLASSES; k++) set(1, k, -1.0f - k);
        for (int k = 0; k < CLASSES; k++) set(2, k, Float.NaN);
        set(3, 2, Float.POSITIVE_INFINITY);
        set(4, 1, 0.75f); set(4, 2, 0.75f);
        set(5, 3, 0.75f); set(5, 2, 0.25f);
        set(6, 0, Float.NaN); set(6, 1, 0.2f);
        set(7, 1, Float.POSITIVE_INFINITY); set(7, 3, Float.POSITIVE_INFINITY);
        set(8, 0, 0.5f); set(8, 1, 0.8f); set(8, 2, 0.8f); set(8, 3, 0.1f);
        BEST[0][0][0] = 1; BEST[0][0][1] = 2; BEST[0][0][2] = 3; BEST[0][0][3] = 0;
    }

    private static byte[][][] generateClassification(boolean useApprox) {
        float max, score;
        int b = 0;
        byte[][][] img = new byte[NX][NY][NZ];
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            if (useApprox) {
                img[x][y][z] = TEMPLATE[BEST[x][y][z]];
            } else {
                // Source excerpt: max starts at 0, b at -1, strict > preserves first ties.
                max = 0; b = -1;
                for (int k = 0; k < CLASSES; k++) {
                    score = MEMS[x][y][z][k];
                    if (score > max) {
                        b = k;
                        max = score;
                    }
                }
                if (b > -1) img[x][y][z] = TEMPLATE[b];
                else img[x][y][z] = 0;
            }
        }
        return img;
    }

    private static String jsonArray(byte[][][] data) {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int z = 0; z < NZ; z++) for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
            if (!first) out.append(',');
            first = false;
            out.append(data[x][y][z]);
        }
        return out.append(']').toString();
    }

    private static String byteArray(byte[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        return out.append(']').toString();
    }

    private static String membershipArray() {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int z = 0; z < NZ; z++) for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++)
            for (int k = 0; k < CLASSES; k++) {
                if (!first) out.append(',');
                first = false;
                float v = MEMS[x][y][z][k];
                if (Float.isNaN(v)) out.append("\"NaN\"");
                else if (v == Float.POSITIVE_INFINITY) out.append("\"Infinity\"");
                else if (v == Float.NEGATIVE_INFINITY) out.append("\"-Infinity\"");
                else out.append(Float.toString(v));
            }
        return out.append(']').toString();
    }

    private static String bestArray() {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int z = 0; z < NZ; z++) for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
            if (!first) out.append(',');
            first = false;
            out.append(BEST[x][y][z]);
        }
        return out.append(']').toString();
    }

    public static void main(String[] args) {
        initialize();
        System.out.print("{\"source_sha256\":\"" + SOURCE_SHA + "\",\"shape\":[" + NX + "," + NY + "," + NZ
            + "],\"template_labels\":" + byteArray(TEMPLATE)
            + ",\"memberships_x_fastest_class_inner\":" + membershipArray()
            + ",\"best_indices_x_fastest\":" + bestArray()
            + ",\"exact_x_fastest\":" + jsonArray(generateClassification(false))
            + ",\"approx_x_fastest\":" + jsonArray(generateClassification(true)) + "}");
    }
}
