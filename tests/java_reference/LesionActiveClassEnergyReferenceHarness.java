/*
 * Selected VandyCRUISE LesionToadSegmentation.computeMemberships() data-term
 * dispatch, source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * Descriptor class order is cruise-atlas-12obj-lesiontoads2012.txt. This
 * evaluates the active OBJECT, LESION, MASK, and fallback branches only;
 * optimized-CSF's old PowerTable implementation is inside a source comment.
 */
public final class LesionActiveClassEnergyReferenceHarness {
    static final String[] TYPES = {
        "mask", "optimized_csf", "optimized_csf", "obj", "obj", "obj",
        "obj", "obj", "obj", "obj", "lesion", "obj", "obj"
    };
    static float[] image = {0.0f, 2.0f};
    static float[] weights = {0.5f, 1.5f};
    static float outlier = 3.25f;
    static float lesionCentroid(int c) { return c == 0 ? 0.25f : 0.75f; }
    static float blackHoleCentroid(int c) { return c == 0 ? 0.5f : 1.25f; }
    static float centroid(int k, int c) { return (k * 2 + c) / 8.0f; }
    static long bits(float value) { return Integer.toUnsignedLong(Float.floatToIntBits(value)); }

    public static void main(String[] args) {
        long[] energies = new long[TYPES.length];
        float lesionEnergy = 0.0f;
        float blackHoleEnergy = 0.0f;
        for (int k = 0; k < TYPES.length; k++) {
            String type = TYPES[k];
            float num = 0.0f, numL = 0.0f, numB = 0.0f;
            if (type.equals("obj") || type.equals("lesion") || type.equals("blackhole")) {
                for (int c = 0; c < image.length; c++) {
                    float d = image[c] - centroid(k, c);
                    num += weights[c] * d * d;
                    if (type.equals("lesion") || type.equals("blackhole")) {
                        float dl = image[c] - lesionCentroid(c);
                        numL += weights[c] * dl * dl;
                    }
                    if (type.equals("blackhole")) {
                        float db = image[c] - blackHoleCentroid(c);
                        numB += weights[c] * db * db;
                    }
                }
            } else if (type.equals("mask")) {
                for (int c = 0; c < image.length; c++)
                    num += weights[c] * (image[c] == 0.0f ? 0.0f : 1.0f);
            } else {
                // Selected active source fallback for currently commented types.
                num = outlier;
            }
            energies[k] = bits(num);
            if (type.equals("lesion") || type.equals("blackhole")) lesionEnergy = numL;
            if (type.equals("blackhole")) blackHoleEnergy = numB;
        }
        System.out.print("{\"class_energy_bits\":[");
        for (int k = 0; k < energies.length; k++) {
            if (k > 0) System.out.print(',');
            System.out.print(energies[k]);
        }
        System.out.println("],\"lesion_energy_bits\":" + bits(lesionEnergy)
            + ",\"black_hole_energy_bits\":" + bits(blackHoleEnergy) + "}");
    }
}
