/** Selected AlgorithmMultipleToads branch after membership and registration. */
public final class ToadsPostRegistrationBranchReferenceHarness {
    private static final String SOURCE_SHA = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final int CURRENT_BITS = 1047831091; // objective from pinned default 11-class reference

    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static void append(StringBuilder out, String name, float previous,
            int iteration, int maximum, float threshold, boolean doGrowing) {
        float current = Float.intBitsToFloat(CURRENT_BITS);
        float ratio = Math.abs(2.0f * (previous-current) / (previous+current));
        boolean stop = false;
        if (ratio < threshold) stop = true;
        if (iteration >= maximum) stop = true;
        boolean runGrowing = (!stop) || doGrowing;
        StringBuilder stages = new StringBuilder("[\"membership_energy\",\"registration_complete\",\"energy_ratio\",\"stop_checks\"");
        if (runGrowing) stages.append(",\"propagateGrowing\"");
        if (stop) stages.append(",\"updateClassificationFromLabels\"");
        else stages.append(",\"computeRelations\",\"computeBoundaryDistances\",\"computeClosestLabels\"");
        stages.append(']');
        if (out.charAt(out.length()-1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"current_bits\":").append(bits(current))
            .append(",\"previous_bits\":").append(bits(previous)).append(",\"ratio_bits\":").append(bits(ratio))
            .append(",\"iteration\":").append(iteration).append(",\"maximum_iterations\":").append(maximum)
            .append(",\"threshold_bits\":").append(bits(threshold)).append(",\"do_growing\":").append(doGrowing)
            .append(",\"stop\":").append(stop).append(",\"run_growing\":").append(runGrowing)
            .append(",\"update_classification\":").append(stop)
            .append(",\"refresh_relations\":").append(!stop).append(",\"growing\":").append(runGrowing)
            .append(",\"stage_order\":").append(stages).append('}');
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"algorithm_source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"post_registration_order\":[\"energyRatio\",\"registerShapesAlreadyComplete\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing\",\"classificationOrRelationRefresh\"],\"cases\":[");
        append(out, "continue-grow-by-source", 0.30f, 1, 5, 0.01f, false);
        append(out, "converged-final-grow", Float.intBitsToFloat(CURRENT_BITS), 2, 5, 0.01f, true);
        append(out, "converged-no-grow", Float.intBitsToFloat(CURRENT_BITS), 2, 5, 0.01f, false);
        append(out, "iteration-limit-grow", 0.30f, 5, 5, 0.01f, true);
        append(out, "iteration-limit-no-grow", 0.30f, 5, 5, 0.01f, false);
        System.out.println(out.append("]}"));
    }
}
