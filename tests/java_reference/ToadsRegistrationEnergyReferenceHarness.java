/*
 * Test-only reference for the exact-membership objective accumulation in
 * VandyCRUISE MultipleToadDeformableAtlas.computeExactRegistrationEnergy.
 * The interpolation result is supplied here as a deterministic synthetic
 * prior field; this excerpt validates source loop order, class selection,
 * positive-prior filtering, float32 arithmetic, normalization, and sampling.
 */
public final class ToadsRegistrationEnergyReferenceHarness {
    private static float[][][][] memberships = new float[5][4][3][4];
    private static float[][][][] priors = new float[5][4][3][4];

    private static float[] evaluate(boolean[] registered, int offset, int subsample) {
        float cost = 0.0f;
        float norm = 0.0f;
        for (int x = offset; x < 5; x += subsample) {
            for (int y = offset; y < 4; y += subsample) {
                for (int z = offset; z < 3; z += subsample) {
                    for (int k = 1; k < 4; k++) {
                        if (registered[k]) {
                            float prior = priors[x][y][z][k];
                            if (prior > 0.0f) {
                                float membership = memberships[x][y][z][k];
                                cost += membership * membership * prior * prior;
                                norm += prior * prior;
                            }
                        }
                    }
                }
            }
        }
        return new float[]{cost, norm, cost / norm};
    }

    private static void print(String name, float[] result, StringBuilder out) {
        if (out.length() > 1) out.append(',');
        out.append('"').append(name).append("\":{")
            .append("\"cost\":").append(Float.toString(result[0])).append(',')
            .append("\"norm\":").append(Float.toString(result[1])).append(',')
            .append("\"energy\":").append(Float.toString(result[2])).append('}');
    }

    public static void main(String[] args) {
        for (int x = 0; x < 5; x++) for (int y = 0; y < 4; y++)
            for (int z = 0; z < 3; z++) for (int k = 0; k < 4; k++) {
                memberships[x][y][z][k] = ((7*x + 3*y + 5*z + 2*k) % 19) / 18.0f;
                priors[x][y][z][k] = ((2*x + 5*y + 3*z + 4*k) % 11) / 10.0f;
            }
        StringBuilder out = new StringBuilder("{");
        print("all_registered", evaluate(new boolean[]{true,true,true,true}, 0, 1), out);
        print("selected_classes_subsampled", evaluate(new boolean[]{false,true,false,true}, 1, 2), out);
        System.out.println(out.append('}'));
    }
}
