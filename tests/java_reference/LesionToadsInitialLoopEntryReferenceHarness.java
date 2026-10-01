import java.util.Locale;

/** Formula excerpt for AlgorithmLesionToads.runAlgorithm() after initial memberships. */
public final class LesionToadsInitialLoopEntryReferenceHarness {
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    // Objective bits emitted by the selected standard-13 volume harness fixture.
    private static final int INITIAL_ENERGY_BITS = 1056474278;

    private static String bits(float value) {
        return Long.toUnsignedString(Integer.toUnsignedLong(Float.floatToRawIntBits(value)));
    }

    private static String oneCase(String name, int iterations, boolean threadStopped, boolean working) {
        float energy = Float.intBitsToFloat(INITIAL_ENERGY_BITS);
        float previousEnergy = energy;
        float energyRatio = 1.0f;
        boolean stop = false;
        int iteration = 0;
        if (iterations == 0) stop = true;
        boolean enter = (!stop) && (!threadStopped) && working;
        if (enter) iteration++;
        return "{\"name\":\"" + name + "\",\"iterations\":" + iterations
            + ",\"thread_stopped\":" + threadStopped + ",\"is_working\":" + working
            + ",\"initial_energy_bits\":" + bits(energy)
            + ",\"previous_energy_bits\":" + bits(previousEnergy)
            + ",\"energy_ratio_bits\":" + bits(energyRatio)
            + ",\"gm_flag_after_initial_memberships\":true"
            + ",\"iteration\":" + iteration + ",\"stop\":" + stop
            + ",\"enter_loop\":" + enter
            + ",\"next_stage\":" + (enter ? "\"propagateThinning\"" : "null") + "}";
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        System.out.print("{\"wrapper_sha256\":\"" + WRAPPER_SHA + "\",\"source_call_order\":["
            + "\"computeMemberships\",\"setGMflag\",\"Eprev=energy\",\"Eratio=1.0f\","
            + "\"stop=false\",\"iter=0\",\"if(iterations==0) stop=true\","
            + "\"while(!stop && !threadStopped && isWorking())\",\"iter++\",\"propagateThinning\"],"
            + "\"cases\":["
            + oneCase("one_iteration_working", 1, false, true) + ","
            + oneCase("zero_iterations", 0, false, true) + ","
            + oneCase("thread_stopped", 1, true, true) + ","
            + oneCase("segmentation_not_working", 1, false, false)
            + "]}");
    }
}
