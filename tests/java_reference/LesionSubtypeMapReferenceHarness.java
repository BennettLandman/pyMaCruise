/**
 * Executable selected-method excerpt for VandyCRUISE
 * LesionToadSegmentation.updateLesionClassifcation().
 *
 * Source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99
 * Wrapper SHA-256: 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb
 *
 * The wrapper calls the update after label transfer when stop is true. The
 * harness preserves the map-clearing else branch and counts qualifying voxels.
 */
public final class LesionSubtypeMapReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final int NX = 7, NY = 1, NZ = 1, CLASSES = 2;
    private static final byte LESION_LABEL = 9;
    private static final byte[] CLASSIFICATION = {9, 9, 9, 9, 9, 8, 9};
    private static final float[] MEMBERSHIP = {
        0.1f, 0.2f,
        0.1f, 0.2f,
        0.1f, 0.2f,
        0.1f, 0.5f,
        0.1f, 0.4f,
        0.1f, 0.1f,
        0.1f, 0.1f
    };
    private static final float[] LESION = {0.4f, 0.3f, 0.5f, 0.5f, 0.3f, 0.8f, 0.1f};
    private static final float[] BLACK_HOLE = {0.3f, 0.5f, 0.5f, 0.5f, 0.39f, 0.7f, 0.9f};

    private static final class Result {
        final short[] map;
        final int count;
        Result(short[] map, int count) { this.map = map; this.count = count; }
    }

    // Faithful selected-method control flow from updateLesionClassifcation().
    private static int updateLesionClassifcation(
            byte[] classification, float[] memberships, float[] lesionMembership,
            float[] blackHoleMembership, boolean hasLesion, short[] lesionClassification) {
        int lesionVolume = 0;
        int lesionClass = 1;
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            int i = (x * NY + y) * NZ + z;
            if (hasLesion)
                if ((classification[i] == LESION_LABEL)
                        && ((lesionMembership[i] > memberships[i * CLASSES + lesionClass])
                            || (blackHoleMembership[i] > memberships[i * CLASSES + lesionClass]))) {
                    if (lesionMembership[i] < blackHoleMembership[i])
                        lesionClassification[i] = (short)2;
                    else
                        lesionClassification[i] = (short)1;
                    lesionVolume++;
                } else
                    lesionClassification[i] = (short)0;
            else
                lesionClassification[i] = (short)0;
        }
        return lesionVolume;
    }

    private static Result run(boolean hasLesion) {
        short[] map = {7, 7, 7, 7, 7, 7, 7};
        int count = updateLesionClassifcation(
            CLASSIFICATION, MEMBERSHIP, LESION, BLACK_HOLE, hasLesion, map);
        return new Result(map, count);
    }

    private static void emit(StringBuilder out, short[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static void emit(StringBuilder out, byte[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static void emit(StringBuilder out, float[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(Float.toString(values[i]));
        }
        out.append(']');
    }

    private static void emitResult(StringBuilder out, String name, Result result) {
        out.append("{\"name\":\"").append(name).append("\",\"lesions_enabled\":")
            .append(name.equals("enabled") ? "true" : "false")
            .append(",\"map\":");
        emit(out, result.map);
        out.append(",\"qualifying_count\":").append(result.count).append('}');
    }

    public static void main(String[] args) {
        Result enabled = run(true);
        Result disabled = run(false);
        StringBuilder out = new StringBuilder("{\"source_method\":\"LesionToadSegmentation.updateLesionClassifcation\"")
            .append(",\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"wrapper_sha256\":\"").append(WRAPPER_SHA)
            .append("\",\"shape\":[").append(NX).append(',').append(NY).append(',').append(NZ)
            .append("],\"lesion_class\":1,\"lesion_label\":").append(LESION_LABEL)
            .append(",\"classification\":");
        emit(out, CLASSIFICATION);
        out.append(",\"template_labels\":[0,").append(LESION_LABEL)
            .append("],\"memberships\":"); emit(out, MEMBERSHIP);
        out.append(",\"lesion_membership\":"); emit(out, LESION);
        out.append(",\"black_hole_membership\":"); emit(out, BLACK_HOLE);
        out.append(",\"initial_map\":[7,7,7,7,7,7,7],\"cases\":[");
        emitResult(out, "enabled", enabled);
        out.append(',');
        emitResult(out, "disabled", disabled);
        System.out.println(out.append("]}").toString());
    }
}
