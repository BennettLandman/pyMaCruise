/** Selected-method excerpt of VandyCRUISE LesionToadSegmentation.updateClassificationFromMemberships. */
public final class LesionClassificationFromMembershipsReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final String WRAPPER_SHA =
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";

    private static byte[] update(byte[] classification, byte[] templateLabels,
            float[][] memberships, float[] lesionMembership, int lesionClass) {
        int classDiff = 0;
        for (int v = 0; v < classification.length; v++) {
            byte prev = classification[v];
            classification[v] = 0;
            int best = 0;
            for (int k = 1; k < templateLabels.length; k++) {
                if (memberships[v][k] > memberships[v][best]) best = k;
            }
            if (lesionMembership[v] > memberships[v][best])
                classification[v] = templateLabels[lesionClass];
            else if (memberships[v][best] > 0)
                classification[v] = templateLabels[best];
            if (prev != classification[v]) classDiff++;
        }
        return new byte[] {(byte) classDiff};
    }

    private static String ints(byte[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        return out.append(']').toString();
    }

    private static String floats(float[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(Float.toString(values[i]));
        }
        return out.append(']').toString();
    }

    private static String memberships(float[][] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(floats(values[i]));
        }
        return out.append(']').toString();
    }

    private static void scenario(StringBuilder out, String name, int[] shape,
            byte[] template, int lesionClass, byte[] previous,
            float[][] memberships, float[] lesionMembership) {
        byte[] result = previous.clone();
        int changed = update(result, template, memberships, lesionMembership, lesionClass)[0];
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[")
            .append(shape[0]).append(',').append(shape[1]).append(',').append(shape[2])
            .append("],\"template_labels\":").append(ints(template))
            .append(",\"lesion_class\":").append(lesionClass)
            .append(",\"previous\":").append(ints(previous))
            .append(",\"memberships\":").append(memberships(memberships))
            .append(",\"lesion_membership\":").append(floats(lesionMembership))
            .append(",\"classification\":").append(ints(result))
            .append(",\"class_diff\":").append(changed).append('}');
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"")
            .append(SOURCE_SHA).append("\",\"wrapper_sha256\":\"")
            .append(WRAPPER_SHA).append("\",\"source_method\":")
            .append("\"LesionToadSegmentation.updateClassificationFromMemberships\",")
            .append("\"scenarios\":[");
        scenario(out, "multiple_classes", new int[] {4, 2, 1},
            new byte[] {5, 37, -7, 99}, 3,
            new byte[] {37, -1, 99, 0, 55, 5, 99, 0},
            new float[][] {
                {0.1f, 0.6f, 0.6f, 0.2f},
                {0f, 0f, 0f, 0f},
                {0.1f, 0.2f, 0.4f, 0.1f},
                {0.1f, 0.3f, 0.7f, 0.2f},
                {0.8f, 0.7f, 0.1f, 0.1f},
                {0.1f, 0.9f, 0.3f, 0.2f},
                {0.1f, 0.2f, 0.5f, 0.5f},
                {0.1f, 0.2f, 0.3f, 0.3f}
            },
            new float[] {0.1f, 0f, 0.5f, 0.7f, 0.2f, 0.91f, 0f, 0f});
        scenario(out, "one_class", new int[] {3, 1, 1},
            new byte[] {-13}, 0,
            new byte[] {-13, 5, 0},
            new float[][] {{0.8f}, {0f}, {0.2f}},
            new float[] {0f, 0.3f, 0.2f});
        System.out.println(out.append("]}"));
    }
}
