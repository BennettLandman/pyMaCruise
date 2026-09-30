/** Selected-source composition of support tiers and their membership factors. */
public final class ToadsRelationFactorCompositionReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int CLASSES = 6;
    private static final byte[][] CLOSEST = {
        {0, 1, 2, 3}, {5, -1, -1, -1}, {-1, 2, 4, -1}, {-1, -1, 3, 3},
        {4, 4, 1, 0}, {-1, -1, -1, -1}, {0, 0, 0, 0}, {2, 3, 2, 3}
    };
    private static final float[] PRIORS = {0.0f, 0.125f, 0.5f, 1.0f, 2.0f};
    private static float[][] shapes() {
        float[][] values = new float[CLOSEST.length][CLASSES];
        for (int v = 0; v < values.length; v++) for (int k = 0; k < CLASSES; k++) {
            if ((v + k) % 4 == 0) values[v][k] = 0.0f;
            else if ((v + k) % 4 == 1) values[v][k] = 1.0e-7f;
            else if ((v + k) % 4 == 2) values[v][k] = -0.25f;
            else values[v][k] = 0.5f;
        }
        return values;
    }
    private static float factor(float shape, float prior, int type) {
        if (type == 2) return 1.0f;
        if (type == 3 && shape == 0.0f) return 0.0f;
        float scaled = 2.0f * prior * shape;
        if (type == 1) {
            float minimum = scaled < 1.0f ? scaled : 1.0f;
            return 0.25f > minimum ? 0.25f : minimum;
        }
        return 1.0f < scaled ? 1.0f : scaled;
    }
    private static void append(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    public static void main(String[] args) {
        float[][] shape = shapes();
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"shape\":[").append(CLOSEST.length).append(",1,1],\"classes\":").append(CLASSES)
            .append(",\"closest_labels\":[");
        boolean first = true;
        for (byte[] row : CLOSEST) for (byte value : row) {
            if (!first) out.append(','); first = false; out.append(value);
        }
        out.append("],\"shape_bits\":["); first = true;
        for (float[] row : shape) for (float value : row) {
            if (!first) out.append(','); first = false; out.append(Float.floatToIntBits(value));
        }
        out.append("],\"priors\":[");
        for (int p = 0; p < PRIORS.length; p++) {
            if (p > 0) out.append(',');
            int[] bits = new int[CLOSEST.length * CLASSES]; int i = 0;
            for (int v = 0; v < CLOSEST.length; v++) for (int k = 0; k < CLASSES; k++) {
                int type = 3;
                if (CLOSEST[v][0] == k || CLOSEST[v][1] == k) type = 2;
                else if (CLOSEST[v][2] == k || CLOSEST[v][3] == k) type = 1;
                bits[i++] = Float.floatToIntBits(factor(shape[v][k], PRIORS[p], type));
            }
            out.append("{\"prior_bits\":").append(Float.floatToIntBits(PRIORS[p])).append(",\"factor_bits\":");
            append(out, bits); out.append('}');
        }
        System.out.println(out.append("]}"));
    }
}
