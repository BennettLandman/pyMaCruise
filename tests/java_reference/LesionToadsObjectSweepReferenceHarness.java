import java.util.Locale;

/** Selected VandyCRUISE computeMemberships OBJECT-only in-place sweep harness. */
public final class LesionToadsObjectSweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx = 5, ny = 4, nz = 6, nc = 3;
        final float[] prior = {0.25f, 0.5f, 0.75f};
        final float smoothing = 0.2f, zero = 1.0e-10f, infinity = 1.0e10f;
        float[][][][] mems = new float[nx][ny][nz][nc];
        float[][][][] data = new float[nx][ny][nz][nc];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            for (int k = 0; k < nc; k++) {
                mems[x][y][z][k] = 0.5f + 0.125f*x + 0.0625f*y + 0.03125f*z + 0.0625f*k;
                data[x][y][z][k] = 0.5f + 0.125f*x + 0.0625f*y + 0.03125f*z + 0.0625f*k;
            }
        }
        int[][] delta = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out = new StringBuilder();
        out.append("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":3,\"smoothing_bits\":").append(bits(smoothing)).append(",\"voxels\":[");
        boolean first = true;
        // Matches computeMemberships: x outermost, then y, then z; immediate writes.
        for (int x = 1; x < nx-1; x++) for (int y = 1; y < ny-1; y++) for (int z = 1; z < nz-1; z++) {
            if (!first) out.append(',');
            first = false;
            float den = 0.0f;
            for (int k = 0; k < nc; k++) {
                float ngb = 0.0f;
                for (int[] face : delta) {
                    for (int m = 0; m < nc; m++) if (m != k) {
                        float value = mems[x+face[0]][y+face[1]][z+face[2]][m];
                        ngb += value*value;
                    }
                }
                float num = data[x][y][z][k];
                num += prior[k]*smoothing*(ngb+0.0f+0.0f);
                if (num > zero) num = 1.0f/num;
                else num = infinity;
                mems[x][y][z][k] = num;
                den += mems[x][y][z][k];
            }
            if (den > 0.0f) for (int k = 0; k < nc; k++) mems[x][y][z][k] = mems[x][y][z][k]/den;
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for (int k = 0; k < nc; k++) {
                if (k > 0) out.append(',');
                out.append(bits(mems[x][y][z][k]));
            }
            out.append("]}");
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
