/*
 * Small executable excerpt of the selected VandyCRUISE
 * MultipleToadSegmentation.generateCruiseOutputs method.
 *
 * Source file SHA-256:
 * 63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884
 * Wrapper file SHA-256:
 * b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d
 *
 * This harness emits tab-separated float32 references for an asymmetric
 * synthetic volume. It intentionally has no dependencies beyond the JDK.
 */
public final class ToadsCruiseOutputsReferenceHarness {
    static final int NX = 3, NY = 4, NZ = 5, N = 14, APPROX = 3;
    static final String[] NAMES = {
        "SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM",
        "CerebralGM", "WhiteMatter", "CerebellarGM", "CerebellumXYZ",
        "Brainstem", "Background", "Sulcalcsf", "Cerebrum-gm",
        "DeepWM", "Cerebellum"
    };

    static float value(int x, int y, int z, int k) {
        return (x * 101 + y * 23 + z * 7 + k * 13 + 1) / 1000.0f;
    }

    static int assigned(int x, int y, int z) {
        return (x + 2 * y + 2 * z) % N;
    }

    static int best(int x, int y, int z, int a) {
        int first = (x + 2 * y + z) % N;
        return (first + new int[] {0, 5, 9}[a]) % N;
    }

    static float[][][][] run(boolean useApprox) {
        float[][][][] memberships = new float[NX][NY][NZ][N];
        byte[][][] segmentation = new byte[NX][NY][NZ];
        int[][][][] best = new int[NX][NY][NZ][APPROX];
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            segmentation[x][y][z] = (byte)assigned(x, y, z);
            for (int k = 0; k < N; k++) memberships[x][y][z][k] = value(x, y, z, k);
            for (int a = 0; a < APPROX; a++) best[x][y][z][a] = best(x, y, z, a);
        }
        float[][][][] out = new float[4][NX][NY][NZ];
        // Excerpt preserves the source's outer atlas-class / inner candidate order.
        for (int k = 0; k < N; k++) {
            if (NAMES[k].equals("SulcalCSF") || NAMES[k].equals("Sulcal-CSF"))
                for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
                    if (useApprox) {
                        for (int a=0;a<APPROX;a++) if (best[x][y][z][a]==k) out[0][x][y][z]=memberships[x][y][z][a];
                    } else out[0][x][y][z]=memberships[x][y][z][k];
                }
        }
        for (int k = 0; k < N; k++) {
            if (NAMES[k].equals("CorticalGM") || NAMES[k].equals("Cerebrum-GM") || NAMES[k].equals("CerebralGM"))
                for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
                    if (useApprox) {
                        for (int a=0;a<APPROX;a++) if (best[x][y][z][a]==k) out[1][x][y][z]=memberships[x][y][z][a];
                    } else out[1][x][y][z]=memberships[x][y][z][k];
                }
        }
        for (int k = 0; k < N; k++) {
            String name = NAMES[k];
            if (!name.equals("SulcalCSF") && !name.equals("Sulcal-CSF")
                    && !name.equals("CorticalGM") && !name.equals("Cerebrum-GM") && !name.equals("CerebralGM")
                    && !name.startsWith("Cerebell") && !name.equals("Brainstem") && !name.equals("Background")) {
                for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
                    if (useApprox) {
                        for (int a=0;a<APPROX;a++) if (best[x][y][z][a]==k) out[2][x][y][z]+=memberships[x][y][z][a];
                    } else out[2][x][y][z]+=memberships[x][y][z][k];
                    if (segmentation[x][y][z]==k) out[3][x][y][z]=1.0f;
                }
            }
        }
        return out;
    }

    static final String[] CHANNELS = {"csf", "gray_matter", "white_matter_fill", "white_matter_mask"};
    public static void main(String[] args) {
        for (boolean approximate : new boolean[] {false, true}) {
            float[][][][] out = run(approximate);
            for (int c=0;c<4;c++) {
                StringBuilder line = new StringBuilder(approximate ? "approx" : "exact").append('\t').append(CHANNELS[c]);
                for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
                    line.append('\t').append(Float.toString(out[c][x][y][z]));
                System.out.println(line);
            }
        }
    }
}
