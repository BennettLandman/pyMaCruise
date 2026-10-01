import java.util.Locale;

/**
 * Raw-bit reference for the selected LesionTOADS energy-ratio update.
 *
 * Wrapper source: VandyCRUISE AlgorithmLesionToads.java, SHA-256
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 * Absolute value: VandyCRUISE Numerics.java, SHA-256
 * 13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb.
 * The method below is a local excerpt of Numerics.abs(float):
 *     if (a > 0) return a;
 *     else return -a;
 *
 * Inputs for the main case are the objective raw bits from the P39ak and
 * P39am selected-Java references. Other cases isolate scalar sign/equality
 * and zero-denominator behavior; this harness does not perform loop control.
 */
public final class LesionEnergyRatioReferenceHarness {
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final String NUMERICS_SHA =
        "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";

    private static float numericsAbs(float value) {
        if (value > 0.0f) return value;
        else return -value;
    }

    private static final class Result {
        final float previous;
        final float current;
        final float difference;
        final float numerator;
        final float denominator;
        final float quotient;
        final float ratio;

        Result(float previous, float current) {
            this.previous = previous;
            this.current = current;
            this.difference = previous - current;
            this.numerator = 2.0f * difference;
            this.denominator = previous + current;
            this.quotient = numerator / denominator;
            this.ratio = numericsAbs(quotient);
        }
    }

    private static String bits(float value) {
        return String.format(Locale.ROOT, "%08x", Float.floatToRawIntBits(value));
    }

    private static void appendCase(StringBuilder out, String name, int previousBits,
                                   int currentBits, boolean first) {
        if (!first) out.append(',');
        Result r = new Result(Float.intBitsToFloat(previousBits), Float.intBitsToFloat(currentBits));
        out.append("{\"case\":\"").append(name).append("\"")
            .append(",\"previous_energy_bits\":\"").append(bits(r.previous)).append("\"")
            .append(",\"current_energy_bits\":\"").append(bits(r.current)).append("\"")
            .append(",\"difference_bits\":\"").append(bits(r.difference)).append("\"")
            .append(",\"numerator_bits\":\"").append(bits(r.numerator)).append("\"")
            .append(",\"denominator_bits\":\"").append(bits(r.denominator)).append("\"")
            .append(",\"quotient_bits\":\"").append(bits(r.quotient)).append("\"")
            .append(",\"energy_ratio_bits\":\"").append(bits(r.ratio)).append("\"}");
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder();
        out.append("{\"wrapper_sha256\":\"").append(WRAPPER_SHA)
            .append("\",\"numerics_sha256\":\"").append(NUMERICS_SHA)
            .append("\",\"prior_objective_reference\":\"lesion_initial_thinning_composition/java_reference.json\"")
            .append(",\"current_objective_reference\":\"lesion_post_thinning_gm_membership_update/java_reference.json\"")
            .append(",\"source_call_order\":[\"Eprev=energy\",\"energy=computeMemberships()\",\"Eratio=Numerics.abs(2.0f*(Eprev-energy)/(Eprev+energy))\"]")
            .append(",\"cases\":[");
        appendCase(out, "p39ak_to_p39am_nan", 0x426d9692, 0x7fc00000, true);
        appendCase(out, "finite_decreasing", 0x424475e0, 0x41d2275e, false);
        appendCase(out, "finite_increasing", 0x41d2275e, 0x424475e0, false);
        appendCase(out, "equal", 0x40400000, 0x40400000, false);
        appendCase(out, "zero_over_zero", 0x00000000, 0x00000000, false);
        appendCase(out, "nonzero_over_zero", 0x3f800000, 0xbf800000, false);
        out.append("]}\n");
        System.out.print(out);
    }
}
