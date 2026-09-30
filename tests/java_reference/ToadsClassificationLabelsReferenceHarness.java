/** Selected-method excerpt of MultipleToadSegmentation.updateClassificationFromLabels. */
public final class ToadsClassificationLabelsReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    private static String array(byte[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i != 0) out.append(',');
            out.append(values[i]);
        }
        return out.append(']').toString();
    }

    private static void emitCase(StringBuilder out, String name, int nx, int ny, int nz,
            byte[] templateLabel, int[] segmentation, byte[] previous) {
        byte[] classification = new byte[segmentation.length];
        int classDiff = 0;
        int index = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            byte prev = previous[index];
            classification[index] = templateLabel[segmentation[index]];
            if (prev != classification[index]) classDiff++;
            index++;
        }
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[")
            .append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"template_labels\":").append(array(templateLabel))
            .append(",\"segmentation\":[");
        for (int i = 0; i < segmentation.length; i++) {
            if (i != 0) out.append(',');
            out.append(segmentation[i]);
        }
        out.append("],\"previous\":").append(array(previous))
            .append(",\"classification\":").append(array(classification))
            .append(",\"class_diff\":").append(classDiff).append('}');
    }

    public static void main(String[] args) {
        byte[] template = new byte[] {0, 2, 4, 17};
        int[] labels = new int[] {
            0, 1, 2, 3, 3, 2, 1, 0,
            1, 1, 2, 2, 3, 3, 0, 0,
            2, 0, 3, 1, 0, 2, 1, 3
        };
        byte[] oldAllZero = new byte[labels.length];
        byte[] oldSame = new byte[] {
            0, 2, 4, 17, 17, 4, 2, 0,
            2, 2, 4, 4, 17, 17, 0, 0,
            4, 0, 17, 2, 0, 4, 2, 17
        };
        oldSame[1] = -1;
        oldSame[10] = 0;
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"")
            .append(SOURCE_SHA).append("\",\"cases\":[");
        emitCase(out, "initialize_from_zero", 3, 2, 4, template, labels, oldAllZero);
        emitCase(out, "partially_changed_with_unchanged_border", 3, 2, 4, template, labels, oldSame);
        System.out.println(out.append("]}"));
    }
}
