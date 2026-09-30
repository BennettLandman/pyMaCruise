/** Source-ordered composition of the selected MultipleToadSegmentation methods. */
public final class ToadsRelationStateRefreshReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static void append(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    private static int[] flattenLabels(byte[][][] a) {
        int[] r = new int[a.length * a[0].length * a[0][0].length]; int p = 0;
        for (byte[][] plane : a) for (byte[] row : plane) for (byte v : row) r[p++] = v;
        return r;
    }
    private static int[] flattenShorts(short[][][] a) {
        int[] r = new int[a.length * a[0].length * a[0][0].length]; int p = 0;
        for (short[][] plane : a) for (short[] row : plane) for (short v : row) r[p++] = v;
        return r;
    }
    private static int[] flattenBytes(byte[][][][] a) {
        int n = a.length * a[0].length * a[0][0].length * a[0][0][0].length;
        int[] r = new int[n]; int p = 0;
        for (byte[][][] plane : a) for (byte[][] row : plane) for (byte[] voxel : row) for (byte v : voxel) r[p++] = v;
        return r;
    }
    private static int[] flattenClosest(byte[][][][] a) {
        int n = a.length * a[0].length * a[0][0].length * 4;
        int[] r = new int[n]; int p = 0;
        for (byte[][][] plane : a) for (byte[][] row : plane) for (byte[] voxel : row) for (byte v : voxel) r[p++] = v;
        return r;
    }
    private static int[] flattenBooleans(boolean[][] a) {
        int[] r = new int[a.length * a[0].length]; int p = 0;
        for (boolean[] row : a) for (boolean v : row) r[p++] = v ? 1 : 0;
        return r;
    }
    public static void main(String[] args) {
        int nx = 9, ny = 8, nz = 7, classes = 6, limit = 2;
        byte[][][] labels = ToadsRelationListReferenceHarness.phantom(nx, ny, nz, classes, 2);
        float[][][][] memberships = ToadsRelationDiffusionReferenceHarness.memberships(nx, ny, nz, classes, false);
        byte[][][][] initial = ToadsBoundaryDistancesReferenceHarness.initial(nx, ny, nz, classes, true);

        // Preserve the wrapper's refresh sequence after a segmentation update.
        ToadsRelationDiffusionReferenceHarness.Outcome relation =
            ToadsRelationDiffusionReferenceHarness.computeRelations(labels, classes, memberships, limit);
        ToadsBoundaryDistancesReferenceHarness.Outcome boundary =
            ToadsBoundaryDistancesReferenceHarness.compute(labels, initial, classes, limit);
        byte[][][][] closest = ToadsClosestLabelsReferenceHarness.closestLabels(labels, limit);
        ToadsRelationListReferenceHarness.Result relationTable =
            ToadsRelationListReferenceHarness.build(labels, classes);

        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"limit\":").append(limit)
            .append(",\"segmentation\":"); append(out, flattenLabels(labels));
        out.append(",\"membership_float_bits\":[");
        boolean first = true;
        for (float[][][] plane : memberships) for (float[][] row : plane) for (float[] voxel : row) for (float v : voxel) {
            if (!first) out.append(','); first = false; out.append(Float.floatToIntBits(v));
        }
        out.append("],\"initial_boundary_distances\":"); append(out, flattenBytes(initial));
        out.append(",\"relation_map\":"); append(out, flattenShorts(relation.map));
        out.append(",\"relation_membership\":"); append(out, flattenBooleans(relationTable.relations));
        out.append(",\"relation_diffusion_iterations\":").append(relation.iterations)
            .append(",\"boundary_diffusion_passes\":").append(boundary.passes)
            .append(",\"boundary_distances\":"); append(out, flattenBytes(boundary.distances));
        out.append(",\"closest_labels\":"); append(out, flattenClosest(closest));
        System.out.println(out.append('}'));
    }
}
