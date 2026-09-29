import java.util.Locale;

/**
 * Focused formula harness for the selected VandyCRUISE
 * LesionToadSegmentation.computeMemberships six-face reads. It preserves the
 * source loop order and face order on coordinate-distinct synthetic arrays;
 * it does not instantiate the full segmentation algorithm.
 */
public final class LesionToadsSixFaceNeighborhoodReferenceHarness {
    private static int bits(float value) {
        return Float.floatToIntBits(value);
    }

    private static void appendBits(StringBuilder out, float[] values) {
        out.append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(bits(values[i]));
        }
        out.append(']');
    }

    private static void appendClassBits(StringBuilder out, float[][] values) {
        out.append('[');
        for (int face = 0; face < values.length; face++) {
            if (face > 0) out.append(',');
            appendBits(out, values[face]);
        }
        out.append(']');
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx = 5, ny = 4, nz = 6, nc = 3;
        float[][][][] mems = new float[nx][ny][nz][nc];
        float[][][] lesion = new float[nx][ny][nz];
        float[][][] blackHole = new float[nx][ny][nz];
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    lesion[x][y][z] = 0.25f + 0.125f * x + 0.0625f * y + 0.03125f * z;
                    blackHole[x][y][z] = 0.5f + 0.0625f * x + 0.03125f * y + 0.015625f * z;
                    for (int c = 0; c < nc; c++) {
                        mems[x][y][z][c] = 1.0f + 0.5f * x + 0.25f * y + 0.125f * z + 0.0625f * c;
                    }
                }
            }
        }

        StringBuilder out = new StringBuilder();
        out.append("{\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz)
           .append("],\"classes\":").append(nc)
           .append(",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\"")
           .append(",\"face_order\":[\"X+\",\"X-\",\"Y+\",\"Y-\",\"Z+\",\"Z-\"],\"voxels\":[");
        boolean first = true;
        // The selected source traverses interior voxels in x/y/z order.
        for (int x = 1; x < nx - 1; x++) {
            for (int y = 1; y < ny - 1; y++) {
                for (int z = 1; z < nz - 1; z++) {
                    if (!first) out.append(',');
                    first = false;
                    int[][] xyz = {{x+1,y,z},{x-1,y,z},{x,y+1,z},{x,y-1,z},{x,y,z+1},{x,y,z-1}};
                    float[][] classFaces = new float[6][nc];
                    float[] lesionFaces = new float[6];
                    float[] blackHoleFaces = new float[6];
                    for (int f = 0; f < 6; f++) {
                        int a = xyz[f][0], b = xyz[f][1], c = xyz[f][2];
                        lesionFaces[f] = lesion[a][b][c];
                        blackHoleFaces[f] = blackHole[a][b][c];
                        for (int k = 0; k < nc; k++) classFaces[f][k] = mems[a][b][c][k];
                    }
                    out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z)
                       .append("],\"class_bits\":");
                    appendClassBits(out, classFaces);
                    out.append(",\"lesion_bits\":");
                    appendBits(out, lesionFaces);
                    out.append(",\"black_hole_bits\":");
                    appendBits(out, blackHoleFaces);
                    out.append('}');
                }
            }
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
