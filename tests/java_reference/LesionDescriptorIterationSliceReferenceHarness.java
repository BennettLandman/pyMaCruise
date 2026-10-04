/*
 * Selected VandyCRUISE LesionToadSegmentation.computeMemberships() excerpt.
 * Source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * Composes active 13-class data-term dispatch with inversion and denominator
 * normalization. All shapes are one, relations are direct, smoothing and
 * prior coefficients are zero, and precomputed statistical state is fixed.
 */
public final class LesionDescriptorIterationSliceReferenceHarness {
    static final String[] TYPES = {
        "mask", "optimized_csf", "optimized_csf", "obj", "obj", "obj",
        "obj", "obj", "obj", "obj", "lesion", "obj", "obj"
    };
    static final float[] IMAGE = {0.0f, 2.0f};
    static final float[] WEIGHT = {0.5f, 1.5f};
    static final float OUTLIER = 3.25f;
    static float centroid(int k, int c) { return (k * 2 + c) / 8.0f; }
    static float lesionCentroid(int c) { return c == 0 ? 0.25f : 0.75f; }
    static float blackHoleCentroid(int c) { return c == 0 ? 0.5f : 1.25f; }
    static long bits(float value) { return Integer.toUnsignedLong(Float.floatToIntBits(value)); }

    public static void main(String[] args) {
        float[] classEnergy = new float[TYPES.length];
        float[] membership = new float[TYPES.length];
        float lesionEnergy = 0.0f;
        float lesionMembership = 0.0f;
        float denominator = 0.0f;
        for (int k = 0; k < TYPES.length; k++) {
            String type = TYPES[k];
            float num = 0.0f, numL = 0.0f;
            if (type.equals("obj") || type.equals("lesion")) {
                for (int c = 0; c < IMAGE.length; c++) {
                    float d = IMAGE[c] - centroid(k, c);
                    num += WEIGHT[c] * d * d;
                    if (type.equals("lesion")) {
                        float dl = IMAGE[c] - lesionCentroid(c);
                        numL += WEIGHT[c] * dl * dl;
                    }
                }
            } else if (type.equals("mask")) {
                for (int c = 0; c < IMAGE.length; c++)
                    num += WEIGHT[c] * (IMAGE[c] == 0.0f ? 0.0f : 1.0f);
            } else {
                num = OUTLIER;
            }
            classEnergy[k] = num;
            if (num > 1.0e-10f) num = 1.0f / num;
            else num = 1.0e10f;
            membership[k] = num;
            if (type.equals("lesion")) {
                lesionEnergy = numL;
                if (numL > 1.0e-10f) numL = 1.0f / numL;
                else numL = 1.0e10f;
                lesionMembership = numL;
                denominator += lesionMembership;
            }
            denominator += membership[k];
        }
        for (int k = 0; k < membership.length; k++) membership[k] /= denominator;
        lesionMembership /= denominator;

        System.out.print("{\"class_energy_bits\":[");
        for (int k = 0; k < classEnergy.length; k++) {
            if (k > 0) System.out.print(',');
            System.out.print(bits(classEnergy[k]));
        }
        System.out.print("],\"class_membership_bits\":[");
        for (int k = 0; k < membership.length; k++) {
            if (k > 0) System.out.print(',');
            System.out.print(bits(membership[k]));
        }
        System.out.println("],\"lesion_energy_bits\":" + bits(lesionEnergy)
            + ",\"lesion_membership_bits\":" + bits(lesionMembership)
            + ",\"denominator_bits\":" + bits(denominator) + "}");
    }
}
