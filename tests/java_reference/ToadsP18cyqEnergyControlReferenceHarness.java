/**
 * Bounded Java-float reference for AlgorithmMultipleToads' energy-ratio and
 * later stop checks. This does not execute the full wrapper or atlas class.
 */
public final class ToadsP18cyqEnergyControlReferenceHarness {
    private static int bits(String text) { return (int)Long.parseLong(text); }
    private static String bool(boolean value) { return value ? "true" : "false"; }
    public static void main(String[] args) {
        if (args.length != 6) throw new IllegalArgumentException(
            "previousBits currentBits thresholdBits iteration maximumIterations doGrowing required");
        float previous = Float.intBitsToFloat(bits(args[0]));
        float current = Float.intBitsToFloat(bits(args[1]));
        float threshold = Float.intBitsToFloat(bits(args[2]));
        int iteration = Integer.parseInt(args[3]);
        int maximumIterations = Integer.parseInt(args[4]);
        boolean doGrowing = Boolean.parseBoolean(args[5]);

        // Exact selected source expression, evaluated before optional atlas
        // registration. Convergence and iteration checks happen afterward.
        float ratio = Math.abs(2.0f * (previous - current) / (previous + current));
        boolean convergenceStop = ratio < threshold;
        boolean iterationLimitStop = iteration >= maximumIterations;
        boolean stop = convergenceStop || iterationLimitStop;
        boolean runGrowing = !stop || doGrowing;
        System.out.print("{\"reference_kind\":\"selected Java float formula\","
            +"\"call_order\":[\"computeExactMemberships\",\"energyRatio\","
            +"\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\"],"
            +"\"ratio_bits\":"+Float.floatToRawIntBits(ratio)
            +",\"convergence_stop\":"+bool(convergenceStop)
            +",\"iteration_limit_stop\":"+bool(iterationLimitStop)
            +",\"stop\":"+bool(stop)+",\"run_growing\":"+bool(runGrowing)+"}");
    }
}
