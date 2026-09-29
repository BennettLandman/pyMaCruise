import java.util.Locale;

/** Selected-method formula harness for streamed six-face smoothing terms. */
public final class LesionToadsVolumeSmoothingReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static void appendUpdates(StringBuilder out, float[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(bits(values[i]));
        }
        out.append(']');
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx = 5, ny = 4, nz = 6, nc = 3;
        final String[] type = {"OBJECT", "LESION", "BLACKHOLE"};
        final float[] prior = {0.25f, 0.5f, 0.75f};
        final float smoothing = 0.2f;
        float[][][][] mems = new float[nx][ny][nz][nc];
        float[][][] lesion = new float[nx][ny][nz];
        float[][][] blackHole = new float[nx][ny][nz];
        float[][][][] classTerm = new float[nx][ny][nz][nc];
        float[][][][] lesionTerm = new float[nx][ny][nz][nc];
        float[][][][] blackHoleTerm = new float[nx][ny][nz][nc];
        float[][][][] outmaskTerm = new float[nx][ny][nz][nc];
        float[][][][] classEnergy = new float[nx][ny][nz][nc];
        float[][][][] lesionEnergy = new float[nx][ny][nz][nc];
        float[][][][] blackHoleEnergy = new float[nx][ny][nz][nc];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            lesion[x][y][z] = 0.25f + 0.125f*x + 0.0625f*y + 0.03125f*z;
            blackHole[x][y][z] = 0.5f + 0.0625f*x + 0.03125f*y + 0.015625f*z;
            for (int k = 0; k < nc; k++) {
                mems[x][y][z][k] = 1.0f + 0.5f*x + 0.25f*y + 0.125f*z + 0.0625f*k;
                classTerm[x][y][z][k] = 0.0625f*(1+x+y+z+k);
                lesionTerm[x][y][z][k] = 0.03125f*(2+x+y+z+k);
                blackHoleTerm[x][y][z][k] = 0.015625f*(3+x+y+z+k);
                outmaskTerm[x][y][z][k] = 0.125f*(1+x+y+z+k);
                classEnergy[x][y][z][k] = 0.25f*(1+x+y+z+k);
                lesionEnergy[x][y][z][k] = 0.125f*(2+x+y+z+k);
                blackHoleEnergy[x][y][z][k] = 0.0625f*(3+x+y+z+k);
            }
        }

        int[][] delta = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out = new StringBuilder();
        out.append("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"smoothing_bits\":").append(bits(smoothing)).append(",\"voxels\":[");
        boolean firstVoxel = true;
        for (int x = 1; x < nx-1; x++) for (int y = 1; y < ny-1; y++) for (int z = 1; z < nz-1; z++) {
            if (!firstVoxel) out.append(',');
            firstVoxel = false;
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"classes\":[");
            for (int k = 0; k < nc; k++) {
                if (k > 0) out.append(',');
                float other = 0, lesionSum = 0, blackHoleSum = 0;
                float[] faces = new float[6];
                float[][] classFaces = new float[6][nc];
                for (int f = 0; f < 6; f++) {
                    int a = x+delta[f][0], b = y+delta[f][1], c = z+delta[f][2];
                    faces[f] = lesion[a][b][c];
                    for (int j = 0; j < nc; j++) {
                        classFaces[f][j] = mems[a][b][c][j];
                        if (j != k) other += classFaces[f][j]*classFaces[f][j];
                    }
                    lesionSum += faces[f]*faces[f];
                    float bh = blackHole[a][b][c];
                    blackHoleSum += bh*bh;
                }
                float current = 0;
                if (type[k].equals("LESION") || type[k].equals("BLACKHOLE")) {
                    for (int f : new int[]{1,0,3,2,5,4}) current += classFaces[f][k]*classFaces[f][k];
                }
                float p = prior[k];
                float all = other + lesionSum;
                all = all + blackHoleSum;
                float full = (p*smoothing)*all;
                float half = ((0.5f*p)*smoothing)*all;
                float cTerm = classTerm[x][y][z][k] + full;
                float cEnergy = classEnergy[x][y][z][k] + half;
                float lTerm = lesionTerm[x][y][z][k];
                float lEnergy = lesionEnergy[x][y][z][k];
                float bTerm = blackHoleTerm[x][y][z][k];
                float bEnergy = blackHoleEnergy[x][y][z][k];
                if (type[k].equals("LESION") || type[k].equals("BLACKHOLE")) {
                    float lSum = other + current;
                    lSum = lSum + blackHoleSum;
                    lTerm += (p*smoothing)*lSum;
                    lEnergy += ((0.5f*p)*smoothing)*lSum;
                }
                if (type[k].equals("BLACKHOLE")) {
                    float bSum = other + current;
                    bSum = bSum + lesionSum;
                    bTerm += (p*smoothing)*bSum;
                    bEnergy += ((0.5f*p)*smoothing)*bSum;
                }
                float[] values = {cTerm,lTerm,bTerm,outmaskTerm[x][y][z][k],cEnergy,lEnergy,bEnergy,other,lesionSum,blackHoleSum,current};
                appendUpdates(out, values);
            }
            out.append("]}");
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
