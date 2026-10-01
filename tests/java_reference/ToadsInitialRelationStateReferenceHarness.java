/** Initial-state reference for AlgorithmMultipleToads relation preparation. */
public final class ToadsInitialRelationStateReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    private static void append(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static int[] flatten(byte[][][] values) {
        int[] result = new int[values.length * values[0].length * values[0][0].length];
        int index = 0;
        for (byte[][] plane : values) for (byte[] row : plane) for (byte value : row)
            result[index++] = value;
        return result;
    }

    private static int[] flatten(short[][][] values) {
        int[] result = new int[values.length * values[0].length * values[0][0].length];
        int index = 0;
        for (short[][] plane : values) for (short[] row : plane) for (short value : row)
            result[index++] = value;
        return result;
    }

    private static int[] flatten(byte[][][][] values) {
        int[] result = new int[values.length * values[0].length * values[0][0].length * values[0][0][0].length];
        int index = 0;
        for (byte[][][] plane : values) for (byte[][] row : plane) for (byte[] voxel : row)
            for (byte value : voxel) result[index++] = value;
        return result;
    }

    private static int[] flatten(boolean[][] values) {
        int[] result = new int[values.length * values[0].length];
        int index = 0;
        for (boolean[] row : values) for (boolean value : row) result[index++] = value ? 1 : 0;
        return result;
    }

    public static void main(String[] args) {
        int nx = 9, ny = 8, nz = 7, classes = 6, limit = 2;
        byte[][][] labels = ToadsRelationListReferenceHarness.phantom(nx, ny, nz, classes, 2);
        float[][][][] memberships = ToadsRelationDiffusionReferenceHarness.memberships(
            nx, ny, nz, classes, false
        );
        byte[][][][] zeroBoundary = ToadsBoundaryDistancesReferenceHarness.initial(
            nx, ny, nz, classes, false
        );

        // Preserve the selected startup call order and constructor-zero boundary state.
        ToadsRelationDiffusionReferenceHarness.Outcome relation =
            ToadsRelationDiffusionReferenceHarness.computeRelations(labels, classes, memberships, limit);
        ToadsBoundaryDistancesReferenceHarness.Outcome boundary =
            ToadsBoundaryDistancesReferenceHarness.compute(labels, zeroBoundary, classes, limit);
        byte[][][][] closest = ToadsClosestLabelsReferenceHarness.closestLabels(labels, limit);
        ToadsRelationListReferenceHarness.Result relationTable =
            ToadsRelationListReferenceHarness.build(labels, classes);

        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"limit\":").append(limit)
            .append(",\"boundary_initialization\":\"zero\",\"segmentation\":");
        append(out, flatten(labels));
        out.append(",\"membership_float_bits\":[");
        boolean first = true;
        for (float[][][] plane : memberships) for (float[][] row : plane) for (float[] voxel : row)
            for (float value : voxel) {
                if (!first) out.append(',');
                first = false;
                out.append(Float.floatToIntBits(value));
            }
        out.append("],\"relation_map\":"); append(out, flatten(relation.map));
        out.append(",\"relation_membership\":"); append(out, flatten(relationTable.relations));
        out.append(",\"relation_diffusion_iterations\":").append(relation.iterations)
            .append(",\"boundary_diffusion_passes\":").append(boundary.passes)
            .append(",\"boundary_distances\":"); append(out, flatten(boundary.distances));
        out.append(",\"closest_labels\":"); append(out, flatten(closest));
        System.out.println(out.append('}'));
    }
}
