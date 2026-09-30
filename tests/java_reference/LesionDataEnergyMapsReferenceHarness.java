/*
 * Selected Java formula excerpt for VandyCRUISE
 * LesionToadSegmentation.computeMemberships(), source SHA-256:
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99
 * Wrapper SHA-256: 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb
 *
 * Scope: unfielded OBJECT/LESION/BLACKHOLE raw squared-residual data terms.
 * The source loops over interior x/y/z voxels, then classes. LESION and
 * BLACKHOLE assign into shared lesion_energy; BLACKHOLE also assigns the
 * shared blackHole_energy. Smoothing, shape priors, inversion, memberships,
 * fields, special CSF handling, and objective accumulation are omitted.
 */
public final class LesionDataEnergyMapsReferenceHarness {
    static final int NX = 5, NY = 6, NZ = 7, NC = 3;
    static final int OBJECT = 1, LESION = 4, BLACKHOLE = 6;

    static float image(int x, int y, int z, int c) {
        return (13 * x - 7 * y + 5 * z + 11 * c - 23) / 8.0f;
    }

    static float centroid(int k, int c) {
        return (2 + 5 * k - 3 * c) / 8.0f;
    }

    static float lesionCentroid(int c) { return (3 + 4 * c) / 8.0f; }
    static float blackHoleCentroid(int c) { return (-5 + 2 * c) / 8.0f; }
    static float weight(int c) { return new float[] {0.5f, 1.25f, 0.75f}[c]; }

    static int bits(float value) { return Float.floatToIntBits(value); }

    static void appendInts(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    static String name(int type) {
        if (type == OBJECT) return "OBJECT";
        if (type == LESION) return "LESION";
        if (type == BLACKHOLE) return "BLACKHOLE";
        throw new IllegalArgumentException("unsupported type " + type);
    }

    static String runCase(String name, int[] types) {
        int voxels = NX * NY * NZ;
        int[] classBits = new int[voxels * types.length];
        int[] lesionByClassBits = new int[voxels * types.length];
        int[] blackHoleByClassBits = new int[voxels * types.length];
        int[] lesionBits = new int[voxels];
        int[] blackHoleBits = new int[voxels];
        int[] centroidBits = new int[types.length * NC];
        for (int k = 0; k < types.length; k++)
            for (int c = 0; c < NC; c++) centroidBits[k * NC + c] = bits(centroid(k, c));

        // Java's selected source uses x outermost, y next, z innermost.
        for (int x = 1; x < NX - 1; x++) for (int y = 1; y < NY - 1; y++) for (int z = 1; z < NZ - 1; z++) {
            int voxel = (x * NY + y) * NZ + z;
            float lesionEnergy = 0.0f, blackHoleEnergy = 0.0f;
            for (int k = 0; k < types.length; k++) {
                float num = 0.0f, numL = 0.0f, numB = 0.0f;
                int type = types[k];
                for (int c = 0; c < NC; c++) {
                    float value = image(x, y, z, c);
                    float d = value - centroid(k, c);
                    num += weight(c) * d * d;
                    if (type == LESION || type == BLACKHOLE) {
                        float dl = value - lesionCentroid(c);
                        numL += weight(c) * dl * dl;
                    }
                    if (type == BLACKHOLE) {
                        float db = value - blackHoleCentroid(c);
                        numB += weight(c) * db * db;
                    }
                }
                classBits[voxel * types.length + k] = bits(num);
                if (type == LESION || type == BLACKHOLE) {
                    lesionByClassBits[voxel * types.length + k] = bits(numL);
                    lesionEnergy = numL;
                }
                if (type == BLACKHOLE) {
                    blackHoleByClassBits[voxel * types.length + k] = bits(numB);
                    blackHoleEnergy = numB;
                }
            }
            lesionBits[voxel] = bits(lesionEnergy);
            blackHoleBits[voxel] = bits(blackHoleEnergy);
        }

        StringBuilder out = new StringBuilder();
        out.append("{\"name\":\"").append(name).append("\",\"class_types\":[");
        for (int k = 0; k < types.length; k++) {
            if (k > 0) out.append(',');
            out.append('"').append(name(types[k])).append('"');
        }
        out.append("],\"class_centroid_bits\":"); appendInts(out, centroidBits);
        out.append(",\"lesion_centroid_bits\":[");
        for (int c = 0; c < NC; c++) { if (c > 0) out.append(','); out.append(bits(lesionCentroid(c))); }
        out.append("],\"black_hole_centroid_bits\":[");
        for (int c = 0; c < NC; c++) { if (c > 0) out.append(','); out.append(bits(blackHoleCentroid(c))); }
        out.append("],\"channel_weight_bits\":[");
        for (int c = 0; c < NC; c++) { if (c > 0) out.append(','); out.append(bits(weight(c))); }
        out.append("],\"class_energy_bits\":"); appendInts(out, classBits);
        out.append(",\"lesion_energy_by_class_bits\":"); appendInts(out, lesionByClassBits);
        out.append(",\"black_hole_energy_by_class_bits\":"); appendInts(out, blackHoleByClassBits);
        out.append(",\"lesion_energy_bits\":"); appendInts(out, lesionBits);
        out.append(",\"black_hole_energy_bits\":"); appendInts(out, blackHoleBits);
        out.append('}');
        return out.toString();
    }

    public static void main(String[] args) {
        String[] cases = {
            runCase("object_lesion_blackhole", new int[] {OBJECT, LESION, BLACKHOLE}),
            runCase("object_blackhole_lesion", new int[] {OBJECT, BLACKHOLE, LESION}),
            runCase("lesion_object_blackhole_lesion", new int[] {LESION, OBJECT, BLACKHOLE, LESION})
        };
        System.out.println("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\","
            + "\"wrapper_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\","
            + "\"shape\":[" + NX + "," + NY + "," + NZ + "],\"channels\":" + NC + ",\"cases\":["
            + String.join(",", cases) + "]}");
    }
}
