/** Selected-method excerpt of MultipleToadSegmentation.computeRelations. */
public final class ToadsRelationDiffusionReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    static final class Outcome {
        short[][][] map;
        int iterations;
    }

    static float[][][][] memberships(int nx, int ny, int nz, int classes, boolean zero) {
        float[][][][] result = new float[nx][ny][nz][classes];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
            for (int k = 0; k < classes; k++) {
                int code = (x * 13 + y * 7 + z * 3 + k * 5) % 23 + 1;
                result[x][y][z][k] = zero ? 0.0f : (float)code / 24.0f;
            }
        return result;
    }

    static Outcome computeRelations(byte[][][] segmentation, int classes,
            float[][][][] memberships, int limit) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        ToadsRelationListReferenceHarness.Result initialized =
            ToadsRelationListReferenceHarness.build(segmentation, classes);
        short[][][] relationMap = initialized.map;
        boolean change = true;
        int iterations = 0;
        for (int t = 1; t < 10 * limit && change; t++) {
            short[][][] tmp = new short[nx][ny][nz];
            change = false;
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
                if (relationMap[x][y][z] < classes) {
                    change = true;
                    float bestMem = 0.0f;
                    tmp[x][y][z] = relationMap[x][y][z];
                    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                        if (i*i+j*j+l*l == 1 && relationMap[x+i][y+j][z+l] >= classes) {
                            float memSum = 0.0f;
                            for (int k = 0; k < classes; k++)
                                if (initialized.relations[relationMap[x+i][y+j][z+l]][k])
                                    memSum += memberships[x+i][y+j][z+l][k];
                            if (memSum > bestMem) {
                                bestMem = memSum;
                                tmp[x][y][z] = relationMap[x+i][y+j][z+l];
                            }
                        }
                    }
                } else tmp[x][y][z] = relationMap[x][y][z];
            }
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++)
                relationMap[x][y][z] = tmp[x][y][z];
            iterations++;
        }
        Outcome outcome = new Outcome(); outcome.map = relationMap; outcome.iterations = iterations; return outcome;
    }

    private static void appendArray(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i != 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }

    private static void emit(StringBuilder out, String name, int nx, int ny, int nz,
            int classes, int mode, int limit, boolean zeroMembership) {
        byte[][][] labels = ToadsRelationListReferenceHarness.phantom(nx, ny, nz, classes, mode);
        float[][][][] mems = memberships(nx, ny, nz, classes, zeroMembership);
        ToadsRelationListReferenceHarness.Result built =
            ToadsRelationListReferenceHarness.build(labels, classes);
        Outcome outcome = computeRelations(labels, classes, mems, limit);
        short[][][] map = outcome.map;
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"limit\":").append(limit)
            .append(",\"zero_membership\":").append(zeroMembership)
            .append(",\"initial_relation_count\":").append(built.relations.length)
            .append(",\"diffusion_iterations\":").append(outcome.iterations).append(",\"segmentation\":[");
        int[] flat = new int[nx * ny * nz]; int p = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) flat[p++] = labels[x][y][z];
        appendArray(out, flat);
        out.append("],\"membership_values\":[");
        p = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) for (int k = 0; k < classes; k++) {
            if (p++ != 0) out.append(','); out.append(Float.floatToIntBits(mems[x][y][z][k]));
        }
        out.append("],\"relation_map\":["); p = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (p++ != 0) out.append(','); out.append(map[x][y][z]);
        }
        out.append("]}");
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        emit(out, "uniform_cap", 7, 7, 7, 4, 0, 1, false);
        emit(out, "planar_zero_limit", 7, 6, 5, 3, 1, 0, false);
        emit(out, "planar_diffusion", 7, 6, 5, 3, 1, 1, false);
        emit(out, "folded_diffusion", 9, 8, 7, 6, 2, 2, false);
        emit(out, "folded_zero_support", 9, 8, 7, 6, 2, 1, true);
        System.out.println(out.append("]}"));
    }
}
