import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Fixed two-pass excerpt of AlgorithmLesionToads.runAlgorithm's active loop.
 * Energy samples are from the P39v/P39w synthetic composition fixture; this
 * harness pins only Java float32 loop controls and source call order. The
 * numeric methods retain their own selected-source harnesses.
 * VandyCRUISE AlgorithmLesionToads.java SHA-256:
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 */
public final class LesionToadsRepeatedIterationReferenceHarness {
    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous - current) / (previous + current));
    }
    private static String quote(String value) { return "\"" + value + "\""; }
    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float[] energies = new float[] {49.18265151977539f, 26.26949119567871f, 25.68988037109375f};
        int nIterMax = 2;
        float maxDistance = 0.0f;
        List<Integer> ratios = new ArrayList<>();
        List<String> reasons = new ArrayList<>();
        List<String> stages = new ArrayList<>();
        for (int iter=1; iter<=Math.abs(nIterMax); iter++) {
            float ePrevious = energies[iter-1];
            float energy = energies[iter];
            float eRatio = ratio(ePrevious, energy);
            ratios.add(Float.floatToRawIntBits(eRatio));
            boolean stop = false;
            if (eRatio < maxDistance) stop = true;
            if (iter >= Math.abs(nIterMax)) stop = true;
            boolean doGrowing = !stop || nIterMax > 0;
            StringBuilder trace = new StringBuilder("[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\"");
            if (doGrowing) trace.append(",\"propagateGrowing\"");
            if (stop) trace.append(",\"updateClassificationFromLabels\",\"updateLesionClassifcation\"]");
            else trace.append(",\"computeRelations\"]");
            stages.add(trace.toString());
            reasons.add(stop ? (eRatio < maxDistance ? "convergence" : "iteration_limit") : "continue");
            if (stop) break;
        }
        System.out.print("{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"energy_bits\":[");
        for (int i=0;i<energies.length;i++) { if (i>0) System.out.print(","); System.out.print(Float.floatToRawIntBits(energies[i])); }
        System.out.print("],\"ratio_bits\":[");
        for (int i=0;i<ratios.size();i++) { if (i>0) System.out.print(","); System.out.print(ratios.get(i)); }
        System.out.print("],\"stop_reasons\":[");
        for (int i=0;i<reasons.size();i++) { if (i>0) System.out.print(","); System.out.print(quote(reasons.get(i))); }
        System.out.print("],\"stage_order\":[");
        for (int i=0;i<stages.size();i++) { if (i>0) System.out.print(","); System.out.print(stages.get(i)); }
        System.out.print("],\"negative_ceiling\":{\"n_iter_max\":-2,\"stop\":true,\"run_growing\":false,\"stop_reason\":\"iteration_limit\"}}");
        System.out.println();
    }
}
