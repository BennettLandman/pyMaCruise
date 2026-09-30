import java.util.Locale;

/**
 * Source excerpt for AlgorithmLesionToads.runAlgorithm() repeated loop
 * control. It records only stop/growth/order decisions; stage numerics are
 * pinned by their individual selected-method Java harnesses.
 * Wrapper SHA-256: 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 */
public final class LesionToadsFixedInputLoopReferenceHarness {
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";

    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous-current) / (previous+current));
    }

    private static String run(String name, float[] energy, float threshold,
                              int startingIteration, int nIterMax) {
        int cap=Math.abs(nIterMax), iteration=startingIteration;
        StringBuilder passes=new StringBuilder("[");
        while(iteration<cap) {
            iteration++;
            float current=energy[iteration-startingIteration];
            float previous=energy[iteration-startingIteration-1];
            float r=ratio(previous,current);
            boolean convergence=r<threshold;
            boolean limit=iteration>=cap;
            boolean stop=convergence||limit;
            boolean growing=!stop||nIterMax>0;
            if(passes.length()>1)passes.append(',');
            passes.append("{\"iteration\":").append(iteration)
                .append(",\"ratio_bits\":").append(Float.floatToRawIntBits(r))
                .append(",\"convergence\":").append(convergence)
                .append(",\"iteration_limit\":").append(limit)
                .append(",\"stop\":").append(stop)
                .append(",\"run_growing\":").append(growing)
                .append(",\"refresh_relations\":").append(!stop)
                .append('}');
            if(stop)break;
        }
        return "{\"name\":\""+name+"\",\"passes\":"+passes.append("]}");
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        System.out.print("{\"wrapper_sha256\":\""+WRAPPER_SHA+
            "\",\"iteration_order\":[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\",\"propagateGrowing\",\"classificationOrRelationRefresh\"],\"cases\":[");
        System.out.print(run("two_continues_then_converges",new float[]{100.0f,80.0f,60.0f,59.0f},0.1f,1,5));
        System.out.print(",");
        System.out.print(run("signed_ceiling_growth_on",new float[]{100.0f,80.0f,60.0f,50.0f},0.0f,1,4));
        System.out.print(",");
        System.out.print(run("signed_ceiling_growth_off",new float[]{100.0f,80.0f,60.0f,50.0f},0.0f,1,-4));
        System.out.println("]}");
    }
}
