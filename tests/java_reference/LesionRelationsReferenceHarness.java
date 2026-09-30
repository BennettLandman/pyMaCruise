import java.util.ArrayList;
import java.util.List;

/** Selected-method excerpts of VandyCRUISE LesionToadSegmentation relations. */
public final class LesionRelationsReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final short EMPTY = -1;
    private static final int MAXLB = 500;

    static final class State {
        short[][][] map;
        boolean[][] list;
        int iterations;
    }

    static byte[][][] phantom(int nx, int ny, int nz, int classes, int mode) {
        byte[][][] labels = new byte[nx][ny][nz];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            int value;
            if (mode == 0) value = (x < nx / 2) ? 0 : 1;
            else if (mode == 1) {
                value = (x * 5 + y * 3 + z * 7 + (x * y + 2 * z) % 5) % classes;
                if (x >= nx / 2 && y > ny / 3 && x + 2 * z > nx + 2) value = (value + 2) % classes;
            } else value = (x + y + z) % classes;
            labels[x][y][z] = (byte)value;
        }
        return labels;
    }

    /** Selected LesionToadSegmentation.buildRelationList excerpt. */
    static State buildRelationList(byte[][][] segmentation, int classes) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        short[][][] relationMap = new short[nx][ny][nz];
        List<short[]> lbs = new ArrayList<>();
        for (short m = 0; m < classes; m++) lbs.add(new short[] {m});
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            relationMap[x][y][z] = segmentation[x][y][z];
        }
        for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
            short[] newlb = new short[MAXLB];
            for (int n = 0; n < MAXLB; n++) newlb[n] = EMPTY;
            int count = 0;
            for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                if (i*i + j*j + l*l == 1 && segmentation[x][y][z] != segmentation[x+i][y+j][z+l]) {
                    short label = segmentation[x+i][y+j][z+l];
                    boolean found = false;
                    for (int b = 0; b < count; b++) if (label == newlb[b]) { found = true; break; }
                    if (!found) newlb[count++] = label;
                }
            }
            if (count > 0) {
                newlb[count++] = segmentation[x][y][z];
                int relationId = -1;
                for (short n = 0; n < lbs.size(); n++) {
                    short[] old = lbs.get(n);
                    boolean same = old.length == count;
                    if (same) for (int b = 0; b < count; b++) {
                        boolean found = false;
                        for (short value : old) if (newlb[b] == value) { found = true; break; }
                        if (!found) { same = false; break; }
                    }
                    if (same) { relationId = n; break; }
                }
                if (relationId < 0) {
                    if (lbs.size() >= MAXLB) throw new IllegalStateException("relation capacity");
                    relationId = lbs.size();
                    short[] relation = new short[count];
                    System.arraycopy(newlb, 0, relation, 0, count);
                    lbs.add(relation);
                }
                relationMap[x][y][z] = (short)relationId;
            } else relationMap[x][y][z] = segmentation[x][y][z];
        }
        boolean[][] relationList = new boolean[lbs.size()][classes];
        for (int n = 0; n < lbs.size(); n++) for (short label : lbs.get(n)) relationList[n][label] = true;
        for (int k = 0; k < classes; k++) for (int n = classes; n < lbs.size(); n++) if (relationList[n][k])
            for (int m = 0; m < classes; m++) if (relationList[n][m]) relationList[k][m] = true;
        State result = new State(); result.map = relationMap; result.list = relationList; return result;
    }

    /** Selected LesionToadSegmentation.computeRelations diffusion excerpt. */
    static State computeRelations(byte[][][] segmentation, int classes, int limit) {
        State state = buildRelationList(segmentation, classes);
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        short[][][] tmp = new short[nx][ny][nz];
        for (int t = 1; t < 5 * limit; t++) {
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
                if (state.map[x][y][z] < classes) {
                    if      (state.map[x-1][y][z] >= classes) tmp[x][y][z] = state.map[x-1][y][z];
                    else if (state.map[x+1][y][z] >= classes) tmp[x][y][z] = state.map[x+1][y][z];
                    else if (state.map[x][y-1][z] >= classes) tmp[x][y][z] = state.map[x][y-1][z];
                    else if (state.map[x][y+1][z] >= classes) tmp[x][y][z] = state.map[x][y+1][z];
                    else if (state.map[x][y][z-1] >= classes) tmp[x][y][z] = state.map[x][y][z-1];
                    else if (state.map[x][y][z+1] >= classes) tmp[x][y][z] = state.map[x][y][z+1];
                    else tmp[x][y][z] = state.map[x][y][z];
                } else tmp[x][y][z] = state.map[x][y][z];
            }
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
                state.map[x][y][z] = tmp[x][y][z];
            }
            state.iterations++;
        }
        return state;
    }

    static void append(StringBuilder out, String name, int nx, int ny, int nz, int classes, int limit, byte[][][] labels) {
        State built = buildRelationList(labels, classes);
        State result = computeRelations(labels, classes, limit);
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"limit\":").append(limit).append(",\"segmentation\":[");
        boolean first = true;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (!first) out.append(','); first = false; out.append(labels[x][y][z]);
        }
        out.append("],\"initial_map\":["); first = true;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (!first) out.append(','); first = false; out.append(built.map[x][y][z]);
        }
        out.append("],\"relation_map\":["); first = true;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (!first) out.append(','); first = false; out.append(result.map[x][y][z]);
        }
        out.append("],\"relation_count\":").append(result.list.length)
            .append(",\"iterations\":").append(result.iterations).append(",\"relations\":["); first = true;
        for (boolean[] row : result.list) for (boolean member : row) {
            if (!first) out.append(','); first = false; out.append(member);
        }
        out.append("]}");
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        append(out, "asymmetric-planar-zero-diffusion", 7, 6, 5, 2, 0, phantom(7, 6, 5, 2, 0));
        append(out, "asymmetric-planar-one-pass", 17, 6, 5, 2, 1, phantom(17, 6, 5, 2, 0));
        append(out, "asymmetric-planar-two-passes", 17, 6, 5, 2, 2, phantom(17, 6, 5, 2, 0));
        append(out, "asymmetric-folded-diffusion", 7, 6, 5, 5, 2, phantom(7, 6, 5, 5, 1));
        append(out, "asymmetric-folded-one-pass", 7, 6, 5, 5, 1, phantom(7, 6, 5, 5, 1));
        append(out, "asymmetric-uniform", 5, 6, 7, 3, 2, phantom(5, 6, 7, 3, 2));
        System.out.println(out.append("]}").toString());
    }
}
