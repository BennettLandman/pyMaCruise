import java.util.Locale;

/**
 * Fixed-input excerpt of the selected VandyCRUISE post-registerShapes branch.
 * Source SHA: AlgorithmLesionToads.java 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 * Numerics.abs(float) is the selected source method, including its NaN sign.
 */
public final class LesionNaNPostRegistrationBranchReferenceHarness {
    private static final class Case {
        final String name;
        final float previous, current, threshold;
        final int iteration, nIterMax;
        Case(String name, float previous, float current, float threshold,
             int iteration, int nIterMax) {
            this.name=name; this.previous=previous; this.current=current;
            this.threshold=threshold; this.iteration=iteration; this.nIterMax=nIterMax;
        }
    }

    private static float ratio(Case c) {
        float value = 2.0f * (c.previous - c.current) / (c.previous + c.current);
        return value > 0.0f ? value : -value; // Numerics.abs(float)
    }

    private static void emit(Case c, boolean first) {
        float eratio = ratio(c);
        boolean convergenceStop = eratio < c.threshold;
        boolean iterationLimitStop = c.iteration >= Math.abs(c.nIterMax);
        boolean stop = convergenceStop || iterationLimitStop;
        boolean runGrowing = !stop || c.nIterMax > 0;
        if (!first) System.out.print(",");
        System.out.printf(Locale.ROOT,
            "{\"name\":\"%s\",\"ratio_bits\":\"%08x\",\"convergence_stop\":%s,\"iteration_limit_stop\":%s,\"stop\":%s,\"run_growing\":%s,\"update_classification\":%s,\"refresh_relations\":%s}",
            c.name, Float.floatToRawIntBits(eratio), convergenceStop,
            iterationLimitStop, stop, runGrowing, stop, !stop);
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float p39akEnergy = Float.intBitsToFloat(0x426d9692);
        float p39amNaNEnergy = Float.intBitsToFloat(0x7fc00000);
        Case[] cases = {
            new Case("nan_before_iteration_limit_negative", p39akEnergy, p39amNaNEnergy, 0.1f, 1, -3),
            new Case("nan_at_iteration_limit_positive", p39akEnergy, p39amNaNEnergy, 0.1f, 3, 3),
            new Case("nan_at_iteration_limit_negative", p39akEnergy, p39amNaNEnergy, 0.1f, 3, -3),
            new Case("finite_threshold_stop", 3.25f, 3.20f, 0.1f, 1, -3),
        };
        System.out.print("{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"order\":[\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing?\",\"classificationOrRelationRefresh\"],\"cases\":[");
        for (int i=0; i<cases.length; i++) emit(cases[i], i==0);
        System.out.println("]}");
    }
}
