import java.util.Locale;

/** Source-formula excerpt of AlgorithmMultipleToads post-membership control. */
public final class ToadsIterationControlReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static void append(StringBuilder out, String name, float previous, float current,
            int iteration, int maximum, float threshold, boolean doGrowing) {
        float ratio = Math.abs(2.0f * (previous-current) / (previous+current));
        boolean stop = false;
        if (ratio < threshold) stop = true;
        if (iteration >= maximum) stop = true;
        boolean runGrowing = (!stop) || doGrowing;
        boolean updateClassification = stop;
        boolean refreshRelations = !stop;
        if (out.charAt(out.length()-1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"previous_bits\":").append(bits(previous))
            .append(",\"current_bits\":").append(bits(current)).append(",\"iteration\":").append(iteration)
            .append(",\"maximum_iterations\":").append(maximum).append(",\"threshold_bits\":").append(bits(threshold))
            .append(",\"do_growing\":").append(doGrowing).append(",\"ratio_bits\":").append(bits(ratio))
            .append(",\"stop\":").append(stop).append(",\"run_growing\":").append(runGrowing)
            .append(",\"update_classification\":").append(updateClassification)
            .append(",\"refresh_relations\":").append(refreshRelations).append('}');
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"cases\":[");
        append(out, "continue", 12.0f, 9.0f, 1, 5, 0.01f, false);
        append(out, "converged-grow", 10.0f, 9.99f, 2, 5, 0.01f, true);
        append(out, "converged-no-grow", 10.0f, 9.99f, 2, 5, 0.01f, false);
        append(out, "limit-grow", 12.0f, 11.0f, 5, 5, 0.01f, true);
        append(out, "limit-no-grow", 12.0f, 11.0f, 5, 5, 0.01f, false);
        System.out.println(out.append("]}"));
    }
}
