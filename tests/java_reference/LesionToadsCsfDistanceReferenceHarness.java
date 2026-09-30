/**
 * Selected-method excerpt for LesionToadSegmentation.computeCSFDistance().
 * The seed selector and byte fastDistance body are copied as a standalone
 * reference; this does not instantiate the full LesionTOADS workflow.
 */
public final class LesionToadsCsfDistanceReferenceHarness {
    private static final String[] NAMES = {
        "Background", "Caudate", "Thalamus", "Putamen", "Subcortical-GM",
        "SubcorticalGM", "BasalGanglia", "Caudate", "caudate", "Thalamus ",
        "Subcortical_GM", "Subcortical-Gm", "Basal-Ganglia", "Other"
    };

    private static final class Case {
        final String name;
        final int nx, ny, nz, cap;
        final int[][][] segmentation;
        Case(String name, int nx, int ny, int nz, int cap) {
            this.name = name; this.nx = nx; this.ny = ny; this.nz = nz; this.cap = cap;
            this.segmentation = new int[nx][ny][nz];
        }
    }

    private static Case[] cases() {
        Case seeded = new Case("all_exact_names_duplicate_and_controls_cap1", 9, 10, 11, 1);
        int[][] seeds = {
            {0, 2, 10, 1}, {1, 3, 8, 2}, {2, 4, 7, 3}, {3, 5, 6, 4},
            {4, 6, 5, 5}, {5, 7, 4, 6}, {8, 9, 0, 7}
        };
        for (int[] p : seeds) seeded.segmentation[p[0]][p[1]][p[2]] = p[3];
        int[][] controls = {
            {1, 8, 10, 8}, {3, 0, 9, 9}, {5, 2, 8, 10},
            {7, 4, 6, 11}, {0, 9, 3, 12}, {8, 0, 5, 13}
        };
        for (int[] p : controls) seeded.segmentation[p[0]][p[1]][p[2]] = p[3];

        Case defaults = new Case("all_exact_names_duplicate_and_controls_default_cap5", 9, 10, 11, 5);
        for (int x = 0; x < seeded.nx; x++) for (int y = 0; y < seeded.ny; y++)
            for (int z = 0; z < seeded.nz; z++) defaults.segmentation[x][y][z] = seeded.segmentation[x][y][z];

        Case empty = new Case("near_match_names_only_empty_mask_default_cap5", 5, 6, 7, 5);
        int[][] near = {
            {0, 0, 0, 8}, {2, 3, 5, 9}, {4, 1, 6, 10},
            {3, 2, 1, 11}, {1, 5, 3, 12}, {4, 5, 0, 13}
        };
        for (int[] p : near) empty.segmentation[p[0]][p[1]][p[2]] = p[3];
        return new Case[] {seeded, defaults, empty};
    }

    private static boolean[][][] computeCsfMask(int[][][] segmentation) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        boolean[][][] mask = new boolean[nx][ny][nz];
        // Exact selected source comparisons; duplicate Caudate indices qualify too.
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            String name = NAMES[segmentation[x][y][z]];
            mask[x][y][z] = name.equals("Caudate") || name.equals("Thalamus")
                || name.equals("Putamen") || name.equals("Subcortical-GM")
                || name.equals("SubcorticalGM") || name.equals("BasalGanglia");
        }
        return mask;
    }

    // Selected signed-byte fastDistance overload, transcribed from the pinned source.
    private static byte[][][] fastDistance(boolean[][][] obj, byte maxDist) {
        int nx = obj.length, ny = obj[0].length, nz = obj[0][0].length;
        byte[][][] dist = new byte[nx][ny][nz];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            dist[x][y][z] = (byte)(obj[x][y][z] ? 0 : -1);
            for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int k = -1; k <= 1; k++) {
                if (i*i + j*j + k*k != 1) continue;
                int a = x+i, b = y+j, c = z+k;
                if (a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz && obj[a][b][c]
                        && dist[x][y][z] == -1) dist[x][y][z] = 1;
            }
        }
        boolean stop = false;
        byte count = 1;
        while (!stop) {
            stop = true;
            for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
                if (dist[x][y][z] != count) continue;
                for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int k = -1; k <= 1; k++) {
                    if (i*i + j*j + k*k != 1) continue;
                    int a = x+i, b = y+j, c = z+k;
                    if (a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz && dist[a][b][c] == -1) {
                        dist[a][b][c] = (byte)(count + 1);
                        stop = false;
                    }
                }
            }
            count = (byte)(count + 1);
            if (count >= maxDist) stop = true;
        }
        return dist;
    }

    private static void flat(StringBuilder out, int[][][] values) {
        out.append('['); boolean first = true;
        for (int[][] plane : values) for (int[] row : plane) for (int value : row) {
            if (!first) out.append(','); first = false; out.append(value);
        }
        out.append(']');
    }
    private static void flat(StringBuilder out, boolean[][][] values) {
        out.append('['); boolean first = true;
        for (boolean[][] plane : values) for (boolean[] row : plane) for (boolean value : row) {
            if (!first) out.append(','); first = false; out.append(value ? 1 : 0);
        }
        out.append(']');
    }
    private static void flat(StringBuilder out, byte[][][] values) {
        out.append('['); boolean first = true;
        for (byte[][] plane : values) for (byte[] row : plane) for (byte value : row) {
            if (!first) out.append(','); first = false; out.append(value);
        }
        out.append(']');
    }

    public static void main(String[] args) {
        Case[] cases = cases();
        StringBuilder out = new StringBuilder("{\"source_method\":\"LesionToadSegmentation.computeCSFDistance\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"default_max_distance\":5,\"wrapper_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"names\":[");
        for (int i = 0; i < NAMES.length; i++) { if (i > 0) out.append(','); out.append('"').append(NAMES[i]).append('"'); }
        out.append("],\"cases\":[");
        for (int i = 0; i < cases.length; i++) {
            Case c = cases[i]; if (i > 0) out.append(',');
            boolean[][][] mask = computeCsfMask(c.segmentation);
            byte[][][] distance = fastDistance(mask, (byte)c.cap);
            out.append("{\"name\":\"").append(c.name).append("\",\"shape\":[")
                .append(c.nx).append(',').append(c.ny).append(',').append(c.nz)
                .append("],\"max_distance\":").append(c.cap).append(",\"segmentation_xyz_flat\":");
            flat(out, c.segmentation); out.append(",\"mask_xyz_flat\":"); flat(out, mask);
            out.append(",\"distance_xyz_flat\":"); flat(out, distance); out.append('}');
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
