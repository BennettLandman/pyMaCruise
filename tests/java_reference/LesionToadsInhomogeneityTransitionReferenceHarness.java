import java.util.Locale;

/**
 * Pinned source-order excerpt from AlgorithmLesionToads.runAlgorithm and
 * InhomogeneityCorrection2.computeExact3DPolynomialField.
 *
 * This harness validates the wrapper's float threshold, one-degree ramp, and
 * per-channel dispatch sequence. The numeric correction kernel is validated
 * separately against its own source reference; this harness does not claim a
 * full workflow comparison.
 */
public final class LesionToadsInhomogeneityTransitionReferenceHarness {
    private static int floatBits(float value) {
        return Float.floatToIntBits(value);
    }

    private static String run(
            String name, float maxDistance, float energyRatio,
            int startingDegree, int polynomialDegree, int channels) {
        final float lowEratio = (float) Math.sqrt(maxDistance);
        int deg = startingDegree;
        if ((deg < polynomialDegree) && (energyRatio < lowEratio)) deg++;
        StringBuilder calls = new StringBuilder();
        for (int c = 0; c < channels; c++) {
            if (c > 0) calls.append(',');
            calls.append(c).append(':').append(deg);
            // InhomogeneityCorrection2.computeExact3DPolynomialField(0) is a no-op.
        }
        return String.format(Locale.ROOT,
            "{\"name\":\"%s\",\"max_distance_bits\":%d,\"energy_ratio_bits\":%d,"
            + "\"low_ratio_bits\":%d,\"starting_degree\":%d,\"degree\":%d,"
            + "\"channels\":%d,\"calls\":\"%s\"}",
            name, floatBits(maxDistance), floatBits(energyRatio),
            floatBits(lowEratio), startingDegree, deg, channels, calls.toString());
    }

    public static void main(String[] args) {
        float threshold = (float) Math.sqrt(0.0001f);
        float below = Math.nextDown(threshold);
        float above = Math.nextUp(threshold);
        System.out.print("{\"source\":\"AlgorithmLesionToads.runAlgorithm degree gate + "
            + "InhomogeneityCorrection2.computeExact3DPolynomialField\","
            + "\"algorithm_source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\","
            + "\"correction_source_sha256\":\"6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8\","
            + "\"cases\":["
            + run("below_threshold", 0.0001f, below, 0, 3, 2) + ","
            + run("at_threshold", 0.0001f, threshold, 1, 3, 2) + ","
            + run("above_threshold", 0.0001f, above, 1, 3, 2) + ","
            + run("advance_one_degree", 0.0001f, 0.004f, 1, 3, 2) + ","
            + run("hold_at_cap", 0.0001f, 0.0f, 3, 3, 2) + ","
            + run("degree_zero_noop_dispatch", 0.0001f, 0.02f, 0, 3, 2) + ","
            + run("nan_ratio_holds", 0.0001f, Float.intBitsToFloat(0x7fc00000), 1, 3, 2)
            + "]}");
    }
}
