/**
 * Minimal source-expression reference for MultipleToadAtlas's exact scaled
 * registration line search and its caller's `diff > minEdiff` loop condition.
 * The selected Java method evaluates a NaN trial energy, rejects it because
 * `NaN > E0` is false, and returns a NaN relative change after its rejection
 * limit. The caller then exits because `NaN > minEdiff` is false.
 */
public final class ToadsNonfiniteRegistrationTrialReferenceHarness {
    private static float energy(float[] transform, float currentEnergy) {
        for (float value : transform) {
            if (!Float.isFinite(value)) {
                // The source interpolation returns no positive prior support
                // at non-finite mapped coordinates, so cost/norm is 0/0.
                return 0.0f / 0.0f;
            }
        }
        return currentEnergy;
    }

    public static void main(String[] args) {
        float[] transform = new float[6];
        float[] gradient = new float[6];
        float[] hessian = new float[6];
        float initialEnergy = 0.25f;
        float value = initialEnergy;
        float lambda = 1.0f;
        boolean accepted = false;
        int rejected = 0;
        boolean stop = false;
        while (!stop) {
            float[] trial = new float[6];
            for (int n = 0; n < trial.length; n++) {
                trial[n] = transform[n] + lambda / hessian[n] * gradient[n];
            }
            value = energy(trial, initialEnergy);
            if (value > initialEnergy) {
                transform = trial;
                lambda = lambda * 1.5f;
                accepted = true;
                stop = true;
            } else {
                lambda = lambda / 1.5f;
                rejected++;
                if (rejected > 10) stop = true;
            }
        }
        float diff = (value - initialEnergy) / value;
        boolean outerContinues = diff > 1e-6f;
        System.out.printf(
            "{\"accepted\":%s,\"rejected_trials\":%d,\"energy_nan\":%s,"
                + "\"relative_change_nan\":%s,\"outer_continues\":%s,"
                + "\"transform_unchanged\":%s}%n",
            accepted,
            rejected,
            Float.isNaN(value),
            Float.isNaN(diff),
            outerContinues,
            !accepted
        );
    }
}
