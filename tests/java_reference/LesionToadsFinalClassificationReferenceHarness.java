/** Selected-method excerpts of the LesionTOADS final label/subtype updates. */
public final class LesionToadsFinalClassificationReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final int NX = 3, NY = 2, NZ = 2, CLASSES = 3;
    private static final byte[] TEMPLATE_LABELS = {0, 15, 22};
    private static final byte[] SEGMENTATION = {1,1,1,0,2,1,1,1,1,1,1,0};
    private static final byte[] PREVIOUS = {0,15,22,0,0,15,0,22,15,22,0,0};
    private static final float[] MEMBERSHIPS = {
        .2f,.2f,.1f, .1f,.4f,.2f, .1f,.5f,.2f, .7f,.1f,.1f,
        .1f,.1f,.6f, .1f,.29f,.1f, .1f,.61f,.1f, .1f,.7f,.1f,
        .1f,.2f,.1f, .1f,.9f,.1f, .1f,.3f,.1f, .8f,.1f,.1f
    };
    private static final float[] LESION = {.4f,.1f,.5f,.9f,.1f,.3f,.61f,.69f,.2f,.89f,.4f,.1f};
    private static final float[] BLACK_HOLE = {.1f,.5f,.5f,.1f,.9f,.3f,.59f,.71f,.9f,.88f,.4f,.9f};

    private static final class Result {
        final byte[] classification;
        final short[] lesionClassification;
        final int classDiff;
        final int lesionVolume;
        Result(byte[] c, short[] l, int d, int v) {
            classification = c; lesionClassification = l; classDiff = d; lesionVolume = v;
        }
    }

    // Mirrors updateClassificationFromLabels followed by the source's
    // updateLesionClassifcation method body, including strict comparisons.
    private static Result execute(boolean hasLesion) {
        byte[] classification = new byte[SEGMENTATION.length];
        int classDiff = 0, lesionVolume = 0, index = 0;
        short[] lesionClassification = new short[SEGMENTATION.length];
        for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) {
            byte prev = PREVIOUS[index];
            classification[index] = TEMPLATE_LABELS[SEGMENTATION[index]];
            if (prev != classification[index]) classDiff++;
            int base = index * CLASSES;
            if (hasLesion && classification[index] == TEMPLATE_LABELS[1]
                    && (LESION[index] > MEMBERSHIPS[base+1]
                        || BLACK_HOLE[index] > MEMBERSHIPS[base+1])) {
                lesionClassification[index] = LESION[index] < BLACK_HOLE[index]
                    ? (short)2 : (short)1;
                lesionVolume++;
            } else {
                lesionClassification[index] = 0;
            }
            index++;
        }
        return new Result(classification, lesionClassification, classDiff, lesionVolume);
    }

    private static void emitArray(StringBuilder out, short[] values) {
        out.append('[');
        for (int i=0; i<values.length; i++) {
            if (i>0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static void emitArray(StringBuilder out, byte[] values) {
        out.append('[');
        for (int i=0; i<values.length; i++) {
            if (i>0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static void emitArray(StringBuilder out, float[] values) {
        out.append('[');
        for (int i=0; i<values.length; i++) {
            if (i>0) out.append(',');
            out.append(Float.toString(values[i]));
        }
        out.append(']');
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"")
            .append(SOURCE_SHA).append("\",\"shape\":[3,2,2],\"template_labels\":[0,15,22],\"inputs\":{")
            .append("\"segmentation\":");
        emitArray(out, SEGMENTATION);
        out.append(",\"previous_classification\":");
        emitArray(out, PREVIOUS);
        out.append(",\"memberships\":");
        emitArray(out, MEMBERSHIPS);
        out.append(",\"lesion_membership\":");
        emitArray(out, LESION);
        out.append(",\"black_hole_membership\":");
        emitArray(out, BLACK_HOLE);
        out.append("},\"cases\":[");
        for (int mode=0; mode<2; mode++) {
            if (mode>0) out.append(',');
            Result result = execute(mode==0);
            out.append("{\"has_lesion\":").append(mode==0)
                .append(",\"classification\":");
            emitArray(out, result.classification);
            out.append(",\"lesion_classification\":");
            emitArray(out, result.lesionClassification);
            out.append(",\"class_diff\":").append(result.classDiff)
                .append(",\"lesion_volume\":").append(result.lesionVolume).append('}');
        }
        System.out.println(out.append("]}").toString());
    }
}
