/** Selected AlgorithmMultipleToads post-membership loop excerpt for P18cyt composition. */
public final class ToadsP18cytRegisteredIterationReferenceHarness {
    private static final String SOURCE_SHA = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";

    public static void main(String[] args) {
        if (args.length != 7) throw new IllegalArgumentException(
            "previousBits currentBits thresholdBits iteration maximum doGrowing alignShapes required");
        float previous = Float.intBitsToFloat(Integer.parseInt(args[0]));
        float current = Float.intBitsToFloat(Integer.parseInt(args[1]));
        float threshold = Float.intBitsToFloat(Integer.parseInt(args[2]));
        int iteration = Integer.parseInt(args[3]);
        int maximum = Integer.parseInt(args[4]);
        boolean doGrowing = Boolean.parseBoolean(args[5]);
        boolean alignShapes = Boolean.parseBoolean(args[6]);

        // Preserve the selected wrapper expression and its executable order.
        float ratio = Math.abs(2.0f * (previous - current) / (previous + current));
        boolean stop = ratio < threshold;
        if (iteration >= maximum) stop = true;
        boolean runGrowing = (!stop) || doGrowing;
        StringBuilder order = new StringBuilder("[\"computeExactMemberships\",\"energyRatio\"");
        if (alignShapes) order.append(",\"atlas.registerShapes\"");
        order.append(",\"convergenceCheck\",\"iterationLimitCheck\"");
        if (runGrowing) order.append(",\"propagateGrowing\"");
        if (stop) order.append(",\"updateClassificationFromLabels\"");
        else order.append(",\"computeRelations\",\"computeBoundaryDistances\",\"computeClosestLabels\"");
        order.append(']');
        System.out.println("{\"source\":\"AlgorithmMultipleToads selected loop after computeMemberships\","
            +"\"source_sha256\":\""+SOURCE_SHA+"\",\"ratio_bits\":"+Float.floatToIntBits(ratio)
            +",\"stop\":"+stop+",\"run_growing\":"+runGrowing+",\"order\":"+order+"}");
    }
}
