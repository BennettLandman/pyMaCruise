/**
 * Pinned selected-method excerpt of VandyCRUISE
 * LesionToadSegmentation.generateClassification().
 *
 * The source excerpt is intentionally kept small and uses fixed synthetic
 * memberships to expose strict tie, nonpositive, NaN, and infinity behavior.
 */
public final class LesionToadsHardClassificationReferenceHarness {
    private static byte[][][] generateClassification(
            float[][][][] mems, float[][][] lesionMembership,
            byte[] templateLabel, int classes, int lesionClass) {
        int nx = mems.length, ny = mems[0].length, nz = mems[0][0].length;
        byte[][][] img = new byte[nx][ny][nz];
        float dist, max, count;
        int best = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            max = 0; best = -1;
            for (int k = 0; k < classes; k++) {
                if (mems[x][y][z][k] > max) {
                    best = k;
                    max = mems[x][y][z][k];
                }
            }
            if (lesionMembership[x][y][z] > max)
                img[x][y][z] = templateLabel[lesionClass];
            else if (best > -1) img[x][y][z] = templateLabel[best];
            else img[x][y][z] = 0;
        }
        return img;
    }

    public static void main(String[] args) {
        float[][][][] mems = new float[10][1][1][3];
        float[][][] lesion = new float[10][1][1];
        float[][] values = {
            {0.2f, 0.4f, 0.3f}, {0.5f, 0.5f, 0.1f}, {0, 0, 0},
            {-1, -2, -3}, {Float.NaN, 0.1f, 0},
            {Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, 0},
            {0.8f, 0.2f, 0}, {0, 0, 0}, {0, 0.4f, 0.4f}, {0, 0, 0}
        };
        float[] lesionValues = {0.5f, 0.5f, 0, -0.5f, 0.2f,
            Float.POSITIVE_INFINITY, 0.8f, 0.01f, 0.1f, Float.NaN};
        for (int x = 0; x < 10; x++) {
            for (int k = 0; k < 3; k++) mems[x][0][0][k] = values[x][k];
            lesion[x][0][0] = lesionValues[x];
        }
        byte[] labels = {10, 20, 30, 99};
        byte[][][] out = generateClassification(mems, lesion, labels, 3, 3);
        System.out.print("{\"source_method\":\"LesionToadSegmentation.generateClassification\",");
        System.out.print("\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",");
        System.out.print("\"template_labels\":[10,20,30,99],\"lesion_class\":3,\"labels\":[");
        for (int x = 0; x < 10; x++) {
            if (x > 0) System.out.print(",");
            System.out.print(out[x][0][0]);
        }
        System.out.println("]}");
    }
}
