
/**
 * Selected-method excerpt harness for LesionToadSegmentation's exact
 * computeBstemDistance name comparisons and byte fastDistance output.
 * This avoids constructing the full MIPAV/JIST segmentation application.
 */
public final class LesionToadsBrainstemDistanceReferenceHarness {
    private static final String[] NAMES = {
        "Background", "Brainstem", "InterVentricular-WM", "FornixWM",
        "brainstem", "InterVentricularWM", "Fornix-WM", "Other"
    };

    private static final class Case {
        final String name;
        final int nx, ny, nz, cap;
        int[][][] segmentation;
        Case(String name, int nx, int ny, int nz, int cap) {
            this.name = name; this.nx = nx; this.ny = ny; this.nz = nz; this.cap = cap;
            this.segmentation = new int[nx][ny][nz];
        }
    }

    private static Case[] cases() {
        Case seeds = new Case("asymmetric_seed_labels_cap4", 7, 8, 9, 4);
        seeds.segmentation[0][2][7] = 1; // border Brainstem
        seeds.segmentation[4][1][3] = 2; // interior InterVentricular-WM
        seeds.segmentation[5][5][2] = 3; // interior FornixWM
        seeds.segmentation[2][6][1] = 4; // exact-name negative control
        seeds.segmentation[1][4][8] = 5;
        seeds.segmentation[6][0][3] = 6;

        Case capOne = new Case("asymmetric_seed_labels_cap1", 7, 8, 9, 1);
        capOne.segmentation = seeds.segmentation;

        Case empty = new Case("alias_only_empty_mask_cap3", 5, 6, 7, 3);
        empty.segmentation[0][0][0] = 4;
        empty.segmentation[2][3][5] = 5;
        empty.segmentation[4][1][6] = 6;
        return new Case[] {seeds, capOne, empty};
    }

    private static boolean[][][] computeBrainstemMask(int[][][] segmentation) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        boolean[][][] mask = new boolean[nx][ny][nz];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            String name = NAMES[segmentation[x][y][z]];
            mask[x][y][z] = name.equals("Brainstem")
                    || name.equals("InterVentricular-WM")
                    || name.equals("FornixWM");
        }
        return mask;
    }

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
                    if (a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz
                            && dist[a][b][c] == -1) {
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
        StringBuilder out = new StringBuilder("{\"source_method\":\"LesionToadSegmentation.computeBstemDistance\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"names\":[");
        for (int n = 0; n < NAMES.length; n++) { if (n > 0) out.append(','); out.append('"').append(NAMES[n]).append('"'); }
        out.append("],\"cases\":[");
        for (int ci = 0; ci < cases.length; ci++) {
            Case c = cases[ci]; if (ci > 0) out.append(',');
            boolean[][][] mask = computeBrainstemMask(c.segmentation);
            byte[][][] distance = fastDistance(mask, (byte)c.cap);
            out.append("{\"name\":\"").append(c.name).append("\",\"shape\":[")
               .append(c.nx).append(',').append(c.ny).append(',').append(c.nz)
               .append("],\"max_distance\":").append(c.cap).append(",\"segmentation_xyz_flat\":");
            flat(out, c.segmentation);
            out.append(",\"mask_xyz_flat\":"); flat(out, mask);
            out.append(",\"distance_xyz_flat\":"); flat(out, distance); out.append('}');
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
