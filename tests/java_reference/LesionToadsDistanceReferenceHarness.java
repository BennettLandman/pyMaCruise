/*
 * Test-only Java reference for the selected VandyCRUISE
 * LesionToadSegmentation.fastDistance(boolean[][][], byte, int, int, int)
 * method. The implementation below preserves that method's initialization,
 * six-neighbor seed rule, full-volume layer scans, signed-byte storage, and
 * cap check. It needs no MIPAV/JIST runtime. This validates the support
 * primitive only, not LesionTOADS segmentation equivalence.
 */
import java.util.Locale;

public final class LesionToadsDistanceReferenceHarness {
    private static final class Case {
        final String name;
        final int nx, ny, nz;
        final byte maxDistance;
        final int[][] objectVoxels;

        Case(String name, int nx, int ny, int nz, int maxDistance, int[][] objectVoxels) {
            this.name = name;
            this.nx = nx;
            this.ny = ny;
            this.nz = nz;
            this.maxDistance = (byte) maxDistance;
            this.objectVoxels = objectVoxels;
        }

        boolean[][][] mask() {
            boolean[][][] result = new boolean[nx][ny][nz];
            for (int[] point : objectVoxels) result[point[0]][point[1]][point[2]] = true;
            return result;
        }
    }

    private static final Case[] CASES = {
        new Case("empty_mask", 4, 3, 2, 5, new int[][] {}),
        new Case("full_mask", 3, 2, 4, 5, fullCoordinates(3, 2, 4)),
        new Case("boundary_seed", 5, 4, 3, 4, new int[][] {{0, 0, 0}, {0, 0, 1}}),
        new Case("isolated_cap_1", 7, 5, 5, 1, new int[][] {{3, 2, 2}}),
        new Case("isolated_cap_5", 7, 5, 5, 5, new int[][] {{3, 2, 2}}),
        new Case("isolated_cap_9", 11, 11, 11, 9, new int[][] {{5, 5, 5}}),
        new Case("diagonal_only", 3, 3, 3, 3, new int[][] {{0, 0, 0}}),
    };

    private static int[][] fullCoordinates(int nx, int ny, int nz) {
        int[][] result = new int[nx * ny * nz][3];
        int i = 0;
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++) result[i++] = new int[] {x, y, z};
        return result;
    }

    private byte[][][] fastDistance(boolean[][][] obj, byte maxDist, int nx, int ny, int nz) {
        byte[][][] dist = new byte[nx][ny][nz];

        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (obj[x][y][z])
                dist[x][y][z] = 0;
            else
                dist[x][y][z] = -1;
            for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) if (i * i + j * j + l * l == 1) {
                if ((x + i >= 0) && (x + i < nx) && (y + j >= 0) && (y + j < ny) && (z + l >= 0) && (z + l < nz)
                        && obj[x + i][y + j][z + l]) {
                    if (dist[x][y][z] == -1)
                        dist[x][y][z] = 1;
                }
            }
        }

        boolean stop = false;
        byte count = 1;
        while (!stop) {
            stop = true;
            for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
                if (dist[x][y][z] == count) {
                    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) if (i * i + j * j + l * l == 1) {
                        if ((x + i >= 0) && (x + i < nx) && (y + j >= 0) && (y + j < ny) && (z + l >= 0) && (z + l < nz)
                                && (dist[x + i][y + j][z + l] == -1)) {
                            dist[x + i][y + j][z + l] = (byte) (count + 1);
                            stop = false;
                        }
                    }
                }
            }
            count = (byte) (count + 1);
            if (count >= maxDist) stop = true;
        }
        return dist;
    }

    private static String ints(byte[][][] values) {
        StringBuilder text = new StringBuilder("[");
        boolean first = true;
        for (byte[][] plane : values) for (byte[] row : plane) for (byte value : row) {
            if (!first) text.append(',');
            text.append(value);
            first = false;
        }
        return text.append(']').toString();
    }

    private static String points(int[][] points) {
        StringBuilder text = new StringBuilder("[");
        for (int i = 0; i < points.length; i++) {
            if (i > 0) text.append(',');
            text.append('[').append(points[i][0]).append(',').append(points[i][1]).append(',').append(points[i][2]).append(']');
        }
        return text.append(']').toString();
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        LesionToadsDistanceReferenceHarness harness = new LesionToadsDistanceReferenceHarness();
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.fastDistance(boolean[][][], byte, int, int, int)\",\"cases\":[");
        for (int i = 0; i < CASES.length; i++) {
            Case c = CASES[i];
            if (i > 0) System.out.print(',');
            byte[][][] result = harness.fastDistance(c.mask(), c.maxDistance, c.nx, c.ny, c.nz);
            System.out.print("{\"name\":\"" + c.name + "\",\"shape\":[" + c.nx + "," + c.ny + "," + c.nz + "],\"max_distance\":" + c.maxDistance);
            System.out.print(",\"object_voxels\":" + points(c.objectVoxels));
            System.out.print(",\"distance_xyz_flat\":" + ints(result) + "}");
        }
        System.out.println("]}");
    }
}
