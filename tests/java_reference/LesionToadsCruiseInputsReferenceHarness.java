/**
 * Standalone excerpt harness for the selected
 * LesionToadSegmentation.generateCruiseInputs() method.
 *
 * The output is a fixed Java-float reference for the synthetic RAS NIfTI
 * fixture documented beside it. It intentionally reproduces the source loops,
 * including last-alias-wins and float32 accumulation order.
 */
public final class LesionToadsCruiseInputsReferenceHarness {
    private static final String[] NAMES = {
        "Background", "SulcalCSF", "Sulcal-CSF", "CorticalGM",
        "Cerebrum-GM", "CerebralGM", "Cerebellum-GM", "Cerebellar-WM",
        "Brainstem", "Other", "Lesion", "BlackHole", "WhiteMatter"
    };

    private static float membership(int voxel, int cls) {
        return (float)((voxel + 1) * (cls + 3) % 23) / 23.0f;
    }

    private static float lesion(int voxel) {
        return (float)((voxel + 2) % 11) / 17.0f;
    }

    private static float blackHole(int voxel) {
        return (float)((voxel + 3) % 7) / 29.0f;
    }

    private static boolean excluded(String name) {
        return name.equals("SulcalCSF") || name.equals("Sulcal-CSF")
            || name.equals("CorticalGM") || name.equals("Cerebrum-GM")
            || name.equals("CerebralGM") || name.startsWith("Cerebell")
            || name.equals("Brainstem") || name.equals("Background");
    }

    private static float[][][][] calculate(int nx, int ny, int nz) {
        float[][][][] mems = new float[4][nx][ny][nz];
        int n = 0;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++, n++) {
            for (int k = 0; k < NAMES.length; k++) {
                float value = membership(n, k);
                if (NAMES[k].equals("SulcalCSF") || NAMES[k].equals("Sulcal-CSF"))
                    mems[0][x][y][z] = value;
                if (NAMES[k].equals("CorticalGM") || NAMES[k].equals("Cerebrum-GM")
                        || NAMES[k].equals("CerebralGM"))
                    mems[1][x][y][z] = value;
            }
            mems[2][x][y][z] = lesion(n) + blackHole(n);
            int seg = (n * 7 + 5) % NAMES.length;
            for (int k = 0; k < NAMES.length; k++) {
                if (!excluded(NAMES[k])) {
                    mems[2][x][y][z] += membership(n, k);
                    if (seg == k) mems[3][x][y][z] = 1.0f;
                }
            }
        }
        return mems;
    }

    private static void appendChannel(StringBuilder out, float[][][] channel) {
        out.append('[');
        boolean first = true;
        for (float[][] plane : channel) for (float[] row : plane) for (float value : row) {
            if (!first) out.append(',');
            first = false;
            out.append(Float.floatToRawIntBits(value));
        }
        out.append(']');
    }

    public static void main(String[] args) {
        float[][][][] result = calculate(5, 4, 3);
        StringBuilder out = new StringBuilder("{\"source_method\":\"LesionToadSegmentation.generateCruiseInputs\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape_xyz\":[5,4,3],\"channels\":[\"CSF\",\"GM\",\"WM\",\"hard-WM\"],\"atlas_names\":[");
        for (int i = 0; i < NAMES.length; i++) {
            if (i > 0) out.append(',');
            out.append('"').append(NAMES[i]).append('"');
        }
        out.append("],\"channel_float32_bits_xyz_flat\":[");
        for (int i = 0; i < result.length; i++) {
            if (i > 0) out.append(',');
            appendChannel(out, result[i]);
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
