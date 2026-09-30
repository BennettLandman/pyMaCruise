/** Selected-method excerpt of MultipleToadSegmentation.computeClosestLabels. */
public final class ToadsClosestLabelsReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final byte EMPTY = -1;
    private static final int NX = 9, NY = 8, NZ = 7;

    private static byte[][][] segmentation() {
        byte[][][] labels = new byte[NX][NY][NZ];
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            int value = (x * 7 + y * 5 + z * 3 + (x * y + z) % 3) % 6;
            if (x >= 4 && y >= 3 && x + z > 8) value = (value + 2) % 6;
            labels[x][y][z] = (byte)value;
        }
        return labels;
    }

    static byte[][][][] closestLabels(byte[][][] segmentation, int limit) {
        int nx = segmentation.length, ny = segmentation[0].length, nz = segmentation[0][0].length;
        byte[][][][] closestLabel = new byte[nx][ny][nz][4];
        byte[][][][] tmpLabel = new byte[nx][ny][nz][4];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            closestLabel[x][y][z][0] = segmentation[x][y][z];
            closestLabel[x][y][z][1] = EMPTY;
            closestLabel[x][y][z][2] = EMPTY;
            closestLabel[x][y][z][3] = EMPTY;
            for (int n = 0; n < 4; n++) tmpLabel[x][y][z][n] = closestLabel[x][y][z][n];
        }
        for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
            for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                if (i * i + j * j + l * l == 1 && segmentation[x+i][y+j][z+l] != segmentation[x][y][z]) {
                    byte label = segmentation[x+i][y+j][z+l];
                    if (closestLabel[x][y][z][1] == EMPTY) closestLabel[x][y][z][1] = label;
                    else if (closestLabel[x][y][z][2] == EMPTY && closestLabel[x][y][z][1] != label)
                        closestLabel[x][y][z][2] = label;
                    else if (closestLabel[x][y][z][3] == EMPTY && closestLabel[x][y][z][1] != label
                            && closestLabel[x][y][z][2] != label)
                        closestLabel[x][y][z][3] = label;
                }
            }
            for (int n = 0; n < 4; n++) tmpLabel[x][y][z][n] = closestLabel[x][y][z][n];
        }
        for (int t = 1; t < 5 * limit; t++) {
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++) {
                if (tmpLabel[x][y][z][1] == EMPTY) {
                    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                        if (i*i+j*j+l*l == 1 && closestLabel[x+i][y+j][z+l][1] != EMPTY
                                && tmpLabel[x][y][z][0] != closestLabel[x+i][y+j][z+l][1])
                            tmpLabel[x][y][z][1] = closestLabel[x+i][y+j][z+l][1];
                    }
                } else if (tmpLabel[x][y][z][2] == EMPTY) {
                    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                        if (i*i+j*j+l*l == 1) {
                            byte one = closestLabel[x+i][y+j][z+l][1];
                            byte two = closestLabel[x+i][y+j][z+l][2];
                            if (one != EMPTY && tmpLabel[x][y][z][0] != one && tmpLabel[x][y][z][1] != one)
                                tmpLabel[x][y][z][2] = one;
                            else if (two != EMPTY && tmpLabel[x][y][z][0] != two && tmpLabel[x][y][z][1] != two)
                                tmpLabel[x][y][z][2] = two;
                        }
                    }
                } else if (tmpLabel[x][y][z][3] == EMPTY) {
                    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
                        if (i*i+j*j+l*l == 1) {
                            byte one = closestLabel[x+i][y+j][z+l][1];
                            byte two = closestLabel[x+i][y+j][z+l][2];
                            byte three = closestLabel[x+i][y+j][z+l][3];
                            if (one != EMPTY && tmpLabel[x][y][z][0] != one && tmpLabel[x][y][z][1] != one
                                    && tmpLabel[x][y][z][2] != one)
                                tmpLabel[x][y][z][3] = one;
                            else if (two != EMPTY && tmpLabel[x][y][z][0] != two && tmpLabel[x][y][z][1] != two
                                    && tmpLabel[x][y][z][2] != two)
                                tmpLabel[x][y][z][3] = two;
                            else if (three != EMPTY && tmpLabel[x][y][z][0] != three && tmpLabel[x][y][z][1] != three
                                    && tmpLabel[x][y][z][2] != three)
                                tmpLabel[x][y][z][3] = three;
                        }
                    }
                }
            }
            for (int x = 1; x < nx - 1; x++) for (int y = 1; y < ny - 1; y++) for (int z = 1; z < nz - 1; z++)
                for (int n = 0; n < 4; n++) closestLabel[x][y][z][n] = tmpLabel[x][y][z][n];
        }
        return closestLabel;
    }

    private static void appendIntArray(StringBuilder out, int[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) { if (i != 0) out.append(','); out.append(values[i]); }
        out.append(']');
    }

    private static void emit(StringBuilder out, int limit, byte[][][] labels) {
        byte[][][][] closest = closestLabels(labels, limit);
        if (out.charAt(out.length() - 1) != '[') out.append(',');
        out.append("{\"limit\":").append(limit).append(",\"shape\":[").append(NX).append(',').append(NY).append(',').append(NZ).append("],\"segmentation\":[");
        int[] flatLabels = new int[NX * NY * NZ]; int p = 0;
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) flatLabels[p++] = labels[x][y][z];
        appendIntArray(out, flatLabels);
        out.append("],\"closest_labels\":[");
        p = 0;
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) for (int n = 0; n < 4; n++) {
            if (p++ != 0) out.append(','); out.append(closest[x][y][z][n]);
        }
        out.append("]}");
    }

    public static void main(String[] args) {
        byte[][][] labels = segmentation();
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        emit(out, 0, labels); emit(out, 1, labels); emit(out, 3, labels);
        System.out.println(out.append("]}"));
    }
}
