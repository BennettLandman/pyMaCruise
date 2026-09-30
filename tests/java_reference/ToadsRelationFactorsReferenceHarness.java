/** Selected helper-method excerpts used by MultipleToadSegmentation.computeExactMemberships. */
public final class ToadsRelationFactorsReferenceHarness {
    private static final String SEGMENTATION_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String NUMERICS_SHA =
        "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";
    private static float distanceFactor(byte boundaryDistance, float limit) {
        return 1.0f / (1.0f + boundaryDistance / limit);
    }
    private static float shapeFactor(float shape, float priorFactor) {
        float scaled = 2.0f * priorFactor * shape;
        return 1.0f < scaled ? 1.0f : scaled;
    }
    private static float neighborFactor(float shape, float priorFactor) {
        float scaled = 2.0f * priorFactor * shape;
        float minimum = scaled < 1.0f ? scaled : 1.0f;
        return 0.25f > minimum ? 0.25f : minimum;
    }
    private static void bits(StringBuilder out, float[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(Float.floatToIntBits(values[i]));
        }
        out.append(']');
    }
    public static void main(String[] args) {
        byte[] distances = new byte[] {-128, -64, -1, 0, 1, 2, 3, 4, 32, 127};
        float[] limits = new float[] {0.2f, 1.0f, 3.5f, 15.0f};
        float[] shapes = new float[] {-10.0f, -0.25f, 0.0f, 0.1249f, 0.125f, 0.25f, 0.5f, 1.0f, 3.0f};
        float[] priors = new float[] {0.0f, 0.125f, 0.5f, 1.0f, 2.0f};
        StringBuilder out = new StringBuilder("{\"source_sha256\":{\"MultipleToadSegmentation.java\":\"")
            .append(SEGMENTATION_SHA).append("\",\"Numerics.java\":\"").append(NUMERICS_SHA)
            .append("\"},\"distances\":[");
        for (int n = 0; n < limits.length; n++) {
            if (n > 0) out.append(',');
            float[] result = new float[distances.length];
            int[] input = new int[distances.length];
            for (int i = 0; i < distances.length; i++) {
                input[i] = distances[i]; result[i] = distanceFactor(distances[i], limits[n]);
            }
            out.append("{\"limit_bits\":").append(Float.floatToIntBits(limits[n]))
                .append(",\"input\":[");
            for (int i = 0; i < input.length; i++) { if (i > 0) out.append(','); out.append(input[i]); }
            out.append("],\"factor_bits\":"); bits(out, result); out.append('}');
        }
        out.append("],\"shapes\":[");
        boolean first = true;
        for (float prior : priors) for (float shape : shapes) {
            if (!first) out.append(','); first = false;
            out.append("{\"shape_bits\":").append(Float.floatToIntBits(shape))
                .append(",\"prior_bits\":").append(Float.floatToIntBits(prior))
                .append(",\"shape_factor_bits\":").append(Float.floatToIntBits(shapeFactor(shape, prior)))
                .append(",\"neighbor_factor_bits\":").append(Float.floatToIntBits(neighborFactor(shape, prior))).append('}');
        }
        System.out.println(out.append("]}"));
    }
}
