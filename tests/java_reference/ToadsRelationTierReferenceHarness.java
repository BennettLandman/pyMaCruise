/** Selected exact-membership relation-tier and zero-shape gate excerpt. */
public final class ToadsRelationTierReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int CLASSES = 6;
    private static final byte[][] CLOSEST = {
        {0, 1, 2, 3}, {5, -1, -1, -1}, {-1, 2, 4, -1}, {-1, -1, 3, 3},
        {4, 4, 1, 0}, {-1, -1, -1, -1}, {0, 0, 0, 0}, {2, 3, 2, 3}
    };
    private static int[] ints(int[][] values) {
        int[] result = new int[values.length * values[0].length]; int p = 0;
        for (int[] row : values) for (int value : row) result[p++] = value;
        return result;
    }
    private static void append(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    public static void main(String[] args) {
        int voxels = CLOSEST.length;
        float[][] shapes = new float[voxels][CLASSES];
        int[][] tiers = new int[voxels][CLASSES];
        int[][] eligible = new int[voxels][CLASSES];
        for (int v = 0; v < voxels; v++) for (int k = 0; k < CLASSES; k++) {
            if ((v + k) % 4 == 0) shapes[v][k] = 0.0f;
            else if ((v + k) % 4 == 1) shapes[v][k] = 1.0e-7f;
            else if ((v + k) % 4 == 2) shapes[v][k] = -0.25f;
            else shapes[v][k] = 0.5f;
            int type = 3;
            if (CLOSEST[v][0] == k || CLOSEST[v][1] == k) type = 2;
            else if (CLOSEST[v][2] == k || CLOSEST[v][3] == k) type = 1;
            tiers[v][k] = type;
            eligible[v][k] = (type != 3 || shapes[v][k] != 0.0f) ? 1 : 0;
        }
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"shape\":[").append(voxels).append(",1,1],\"classes\":").append(CLASSES)
            .append(",\"closest_labels\":[");
        boolean first = true;
        for (byte[] row : CLOSEST) for (byte label : row) {
            if (!first) out.append(','); first = false; out.append(label);
        }
        out.append("],\"shape_bits\":["); first = true;
        for (float[] row : shapes) for (float value : row) {
            if (!first) out.append(','); first = false; out.append(Float.floatToIntBits(value));
        }
        out.append("],\"relation_types\":"); append(out, ints(tiers));
        out.append(",\"eligible\":"); append(out, ints(eligible));
        System.out.println(out.append('}'));
    }
}
