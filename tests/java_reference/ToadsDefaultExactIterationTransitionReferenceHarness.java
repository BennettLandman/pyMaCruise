import java.util.Locale;

/** Source-ordered post-membership branch composition from VandyCRUISE. */
public final class ToadsDefaultExactIterationTransitionReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static void append(StringBuilder out, String name, float previous, float current,
            int iteration, int maximum, float threshold, boolean doGrowing) {
        float ratio = Math.abs(2.0f * (previous-current) / (previous+current));
        boolean stop = false;
        if (ratio < threshold) stop = true;
        if (iteration >= maximum) stop = true;
        boolean runGrowing = (!stop) || doGrowing;
        if (out.charAt(out.length()-1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"previous_bits\":")
            .append(bits(previous)).append(",\"current_bits\":").append(bits(current))
            .append(",\"iteration\":").append(iteration)
            .append(",\"maximum_iterations\":").append(maximum)
            .append(",\"threshold_bits\":").append(bits(threshold))
            .append(",\"do_growing\":").append(doGrowing)
            .append(",\"ratio_bits\":").append(bits(ratio))
            .append(",\"stop\":").append(stop)
            .append(",\"run_growing\":").append(runGrowing)
            .append(",\"source_order\":\"")
            .append(stop
                ? (runGrowing ? "membership,control,growing,updateClassificationFromLabels" : "membership,control,updateClassificationFromLabels")
                : "membership,control,growing,computeRelations,computeBoundaryDistances,computeClosestLabels")
            .append("\"}");
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        StringBuilder out = new StringBuilder(
            "{\"algorithm_source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\","
            + "\"segmentation_source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\","
            + "\"fixed_registration\":true,\"cases\":[");
        append(out, "continue-grows-even-when-option-off", 1.0f, 0.5f, 1, 5, 1e-5f, false);
        append(out, "iteration-limit-final-growth", 0.9f, 0.8f, 3, 3, 1e-5f, true);
        append(out, "iteration-limit-no-final-growth", 0.9f, 0.8f, 3, 3, 1e-5f, false);
        System.out.println(out.append("]}").toString());
    }
}
