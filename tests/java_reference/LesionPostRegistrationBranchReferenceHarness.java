import java.util.Locale;

/** Fixed-input Java reference for AlgorithmLesionToads' post-registration branch.
 * Pinned to VandyCRUISE AlgorithmLesionToads.java SHA-256
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 * The active wrapper computes the float32 ratio, checks convergence, then the
 * iteration ceiling, and gates growth using the original signed nIterMax.
 */
public final class LesionPostRegistrationBranchReferenceHarness {
    private static final class Case {
        final String name;
        final float previous, current, threshold;
        final int iteration, nIterMax;
        final boolean narrow;
        Case(String n, float p, float c, float t, int i, int max, boolean nb) {
            name=n; previous=p; current=c; threshold=t; iteration=i; nIterMax=max; narrow=nb;
        }
    }
    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous-current) / (previous+current));
    }
    private static boolean stop(Case c) {
        boolean value = ratio(c.previous,c.current) < c.threshold;
        if (c.iteration >= Math.abs(c.nIterMax)) value=true;
        return value;
    }
    private static boolean grow(Case c) {
        return !c.narrow && (!stop(c) || c.nIterMax > 0);
    }
    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        Case[] cases = new Case[] {
            new Case("continue", 3.25f, 2.5f, 0.1f, 1, 3, false),
            new Case("threshold_stop_final_growth_on", 3.25f, 3.20f, 0.1f, 1, 3, false),
            new Case("threshold_stop_final_growth_off", 3.25f, 3.20f, 0.1f, 1, -3, false),
            new Case("iteration_stop_final_growth_on", 3.25f, 2.5f, 0.1f, 3, 3, false),
            new Case("iteration_stop_final_growth_off", 3.25f, 2.5f, 0.1f, 3, -3, false),
            new Case("narrow_band_suppression", 3.25f, 2.5f, 0.1f, 1, 3, true),
        };
        System.out.print("{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"post_registration_order\":[\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing\",\"classificationOrRelationRefresh\"],\"cases\":[");
        for (int n=0;n<cases.length;n++) {
            Case c=cases[n];
            boolean stop=stop(c);
            boolean grow=grow(c);
            float e=ratio(c.previous,c.current);
            if(n>0) System.out.print(",");
            System.out.printf(Locale.ROOT,
                "{\"name\":\"%s\",\"ratio_bits\":%d,\"stop\":%s,\"run_growing\":%s,\"update_classification\":%s,\"refresh_relations\":%s}",
                c.name, Float.floatToRawIntBits(e), stop, grow, stop, !stop);
        }
        System.out.println("]}");
    }
}
