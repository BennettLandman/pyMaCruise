import java.util.Locale;

/**
 * Selected AlgorithmLesionToads.runAlgorithm() suffix after one continuing
 * P39 iteration: membership/energy, registration, stop controls, then either
 * growth plus final outputs or relation refresh. Wrapper SHA-256:
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 */
public final class LesionContinuingPostMembershipReferenceHarness {
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final String ATLAS_SHA =
        "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e";
    private static final String SEGMENTATION_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";

    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous-current) / (previous+current));
    }

    private static void emit(String name, float previous, float current,
                             float threshold, int iteration, int nIterMax) {
        float energyRatio = ratio(previous, current);
        boolean convergence = energyRatio < threshold;
        boolean iterationLimit = iteration >= Math.abs(nIterMax);
        boolean stop = convergence || iterationLimit;
        boolean runGrowing = !stop || nIterMax > 0;
        System.out.printf(Locale.ROOT,
            "{\"name\":\"%s\",\"ratio_bits\":%d,\"stop\":%s,\"run_growing\":%s,\"final_outputs\":%s,\"refresh_relations\":%s}",
            name, Float.floatToRawIntBits(energyRatio), stop, runGrowing, stop, !stop);
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        System.out.print("{\"wrapper_sha256\":\""+WRAPPER_SHA+
            "\",\"atlas_sha256\":\""+ATLAS_SHA+
            "\",\"segmentation_sha256\":\""+SEGMENTATION_SHA+
            "\",\"order\":[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing\",\"classificationOrRelationRefresh\"],\"controls\":[");
        emit("continue", 3.25f, 2.5f, 0.0f, 2, 4);
        System.out.print(",");
        emit("threshold_stop_growth_off", 3.25f, 3.20f, 0.1f, 2, -4);
        System.out.print(",");
        emit("threshold_stop_growth_on", 3.25f, 3.20f, 0.1f, 2, 4);
        System.out.print(",");
        emit("iteration_stop_growth_off", 3.25f, 2.5f, 0.0f, 4, -4);
        System.out.print(",");
        emit("iteration_stop_growth_on", 3.25f, 2.5f, 0.0f, 4, 4);
        System.out.println("]}");
    }
}
