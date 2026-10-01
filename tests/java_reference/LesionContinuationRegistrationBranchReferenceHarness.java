import java.util.Locale;

/**
 * Fixed-profile wrapper reference for one P39v-result continuation through
 * the selected AlgorithmLesionToads registration and stop/continue decision.
 *
 * The optimizer and post-registration numeric operators are validated in
 * their own Java method harnesses. This excerpt pins their wrapper order and
 * the Java float32 control result for the supplied P39v synthetic profile.
 * Source: VandyCRUISE AlgorithmLesionToads.java SHA-256
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 */
public final class LesionContinuationRegistrationBranchReferenceHarness {
    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous - current) / (previous + current));
    }
    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float previous = 49.18265151977539f;
        float current = 26.26949119567871f;
        float threshold = 0.01f;
        int iteration = 1;
        int nIterMax = 3;
        float energyRatio = ratio(previous, current);
        boolean stop = energyRatio < threshold;
        if (iteration >= Math.abs(nIterMax)) stop = true;
        boolean runGrowing = !stop || nIterMax > 0;
        System.out.printf(Locale.ROOT,
            "{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"previous_energy_bits\":%d,\"membership_energy_bits\":%d,\"ratio_bits\":%d,\"stop\":%s,\"run_growing\":%s,\"refresh_relations\":%s,\"stage_order\":[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing\",\"computeRelations\"]}%n",
            Float.floatToRawIntBits(previous), Float.floatToRawIntBits(current),
            Float.floatToRawIntBits(energyRatio), stop, runGrowing, !stop);
    }
}
