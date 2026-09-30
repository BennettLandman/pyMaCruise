/** Selected-method excerpt of MultipleToadSegmentation.computeBoundaryDistances. */
public final class ToadsBoundaryDistancesReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int[][] FACES = {{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};

    static final class Outcome { byte[][][][] distances; int passes; }

    static Outcome compute(byte[][][] labels, byte[][][][] initial, int classes, int limit) {
        int nx = labels.length, ny = labels[0].length, nz = labels[0][0].length;
        byte[][][][] boundaryDist = new byte[nx][ny][nz][classes];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
            for (int k = 0; k < classes; k++) boundaryDist[x][y][z][k] = initial[x][y][z][k];

        for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++)
            for (int[] face : FACES) {
                int i=face[0], j=face[1], l=face[2];
                if (labels[x][y][z] != labels[x+i][y+j][z+l]) {
                    boundaryDist[x+i][y+j][z+l][labels[x][y][z]] = 1;
                    boundaryDist[x][y][z][labels[x+i][y+j][z+l]] = 1;
                }
            }

        byte[][][] tmp = new byte[nx][ny][nz];
        boolean change = true;
        int passes = 0;
        for (int t = 1; t < 5 * limit && change; t++) {
            for (int k = 0; k < classes; k++) {
                change = false;
                for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
                    tmp[x][y][z] = boundaryDist[x][y][z][k];
                    if (labels[x][y][z] != k && boundaryDist[x][y][z][k] == 0) {
                        for (int[] face : FACES) {
                            int i=face[0], j=face[1], l=face[2];
                            if (boundaryDist[x+i][y+j][z+l][k] > 0)
                                tmp[x][y][z] = (byte)Math.min(tmp[x][y][z], boundaryDist[x+i][y+j][z+l][k] + 1);
                        }
                    }
                }
                for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++)
                    boundaryDist[x][y][z][k] = tmp[x][y][z];
            }
            passes++;
        }
        Outcome outcome = new Outcome(); outcome.distances = boundaryDist; outcome.passes = passes; return outcome;
    }

    static byte[][][][] initial(int nx, int ny, int nz, int classes, boolean prepopulated) {
        byte[][][][] result = new byte[nx][ny][nz][classes];
        if (prepopulated) for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
            for (int k = 0; k < classes; k++) result[x][y][z][k] = (byte)((x + 2*y + 3*z + k) % 4);
        return result;
    }

    private static void appendArray(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i != 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }

    private static void emit(StringBuilder out, String name, int nx, int ny, int nz,
            int classes, int mode, int limit, boolean prepopulated) {
        byte[][][] labels = ToadsRelationListReferenceHarness.phantom(nx, ny, nz, classes, mode);
        byte[][][][] before = initial(nx, ny, nz, classes, prepopulated);
        Outcome result = compute(labels, before, classes, limit);
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"limit\":").append(limit)
            .append(",\"prepopulated\":").append(prepopulated).append(",\"diffusion_passes\":").append(result.passes)
            .append(",\"segmentation\":[");
        int count = nx * ny * nz; int[] flat = new int[count]; int p=0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) flat[p++] = labels[x][y][z];
        appendArray(out, flat); out.append("],\"initial\":["); p=0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) for (int k = 0; k < classes; k++) {
            if (p++ != 0) out.append(','); out.append(before[x][y][z][k]);
        }
        out.append("],\"boundary_distances\":["); p=0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) for (int k = 0; k < classes; k++) {
            if (p++ != 0) out.append(','); out.append(result.distances[x][y][z][k]);
        }
        out.append("]}");
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        emit(out, "planar_no_diffusion", 7, 6, 5, 3, 1, 0, false);
        emit(out, "planar_one_diffusion_pass", 7, 6, 5, 3, 1, 1, false);
        emit(out, "folded_three_limits", 9, 8, 7, 6, 2, 3, false);
        emit(out, "folded_prepopulated", 9, 8, 7, 6, 2, 1, true);
        System.out.println(out.append("]}"));
    }
}
