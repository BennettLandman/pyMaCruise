/*
 * Selected LesionToadSegmentation.computeMemberships() normalization excerpt.
 * Source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * Fixed data energies, direct relations, unit shapes, and zero smoothing/prior
 * isolate ordered class/lesion denominator updates for the bundled descriptor.
 */
public final class LesionDescriptorMembershipReferenceHarness {
    static final String[] TYPES = {
        "mask", "optimized_csf", "optimized_csf", "obj", "obj", "obj",
        "obj", "obj", "obj", "obj", "lesion", "obj", "obj"
    };
    static long bits(float value) { return Integer.toUnsignedLong(Float.floatToIntBits(value)); }

    public static void main(String[] args) {
        float[] memberships = new float[TYPES.length];
        float lesionMembership = 0.0f;
        float denominator = 0.0f;
        float lesionRaw = 2.0f;
        for (int k = 0; k < TYPES.length; k++) {
            float num = 1.0f;
            if (num > 1.0e-10f) num = 1.0f / num;
            else num = 1.0e10f;
            memberships[k] = num; // DIRECT relation, unit shape.
            if (TYPES[k].equals("lesion")) {
                float numL = lesionRaw;
                if (numL > 1.0e-10f) numL = 1.0f / numL;
                else numL = 1.0e10f;
                lesionMembership = numL;
                denominator += lesionMembership;
            }
            denominator += memberships[k];
        }
        if (denominator > 0.0f) {
            for (int k = 0; k < TYPES.length; k++) memberships[k] /= denominator;
            lesionMembership /= denominator;
        }
        System.out.print("{\"denominator_bits\":" + bits(denominator)
            + ",\"class_membership_bits\":[");
        for (int k = 0; k < memberships.length; k++) {
            if (k > 0) System.out.print(',');
            System.out.print(bits(memberships[k]));
        }
        System.out.println("],\"lesion_membership_bits\":" + bits(lesionMembership) + "}");
    }
}
