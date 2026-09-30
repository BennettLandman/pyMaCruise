import java.util.ArrayList;
import java.util.List;

/** Selected-method excerpt of MultipleToadSegmentation.buildRelationList. */
public final class ToadsRelationListReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final short EMPTY = -1;
    private static final int MAXLB = 500;

    static byte[][][] phantom(int nx, int ny, int nz, int classes, int mode) {
        byte[][][] labels = new byte[nx][ny][nz];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            int value;
            if (mode == 0) value = 1;
            else if (mode == 1) value = x < nx / 2 ? 0 : 1;
            else {
                value = (x * 7 + y * 5 + z * 3 + (x * y + z) % 3) % classes;
                if (x >= nx / 2 && y >= ny / 3 && x + z > nx) value = (value + 2) % classes;
            }
            labels[x][y][z] = (byte)value;
        }
        return labels;
    }

    static final class Result {
        short[][][] map;
        boolean[][] relations;
    }

    static Result build(byte[][][] segmentation, int classes) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        short[][][] relationMap = new short[nx][ny][nz];
        List<short[]> lbs = new ArrayList<>();
        for (short m = 0; m < classes; m++) lbs.add(new short[] {m});
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
            relationMap[x][y][z] = segmentation[x][y][z];

        for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
            short[] newlb = new short[MAXLB];
            for (int n = 0; n < MAXLB; n++) newlb[n] = EMPTY;
            int count = 0;
            for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                if (i*i+j*j+l*l == 1 && segmentation[x][y][z] != segmentation[x+i][y+j][z+l]) {
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
        Result result = new Result(); result.map = relationMap; result.relations = relationList; return result;
    }

    private static void append(StringBuilder out, String name, int nx, int ny, int nz, int classes, byte[][][] labels) {
        Result result = build(labels, classes);
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
            .append("],\"classes\":").append(classes).append(",\"segmentation\":[");
        boolean first = true;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (!first) out.append(','); first = false; out.append(labels[x][y][z]);
        }
        out.append("],\"relation_count\":").append(result.relations.length).append(",\"relation_map\":[");
        first = true;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            if (!first) out.append(','); first = false; out.append(result.map[x][y][z]);
        }
        out.append("],\"relations\":["); first = true;
        for (boolean[] row : result.relations) for (boolean member : row) {
            if (!first) out.append(','); first = false; out.append(member);
        }
        out.append("]}");
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        append(out, "uniform", 5, 5, 5, 4, phantom(5, 5, 5, 4, 0));
        append(out, "planar", 7, 6, 5, 3, phantom(7, 6, 5, 3, 1));
        append(out, "folded", 9, 8, 7, 6, phantom(9, 8, 7, 6, 2));
        System.out.println(out.append("]}"));
    }
}
