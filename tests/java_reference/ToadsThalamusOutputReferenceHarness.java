/*
 * Minimal executable excerpt of selected VandyCRUISE
 * MultipleToadSegmentation.generateThalamusOutput().
 *
 * Selected source SHA-256:
 * 63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884
 * Wrapper SHA-256:
 * b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d
 *
 * The input generator creates an asymmetric synthetic array. The two nested
 * loops below preserve the source's atlas-order / candidate-slot overwrite
 * semantics. This excerpt intentionally omits segmentation evolution and
 * MIPAV image writing.
 */
public final class ToadsThalamusOutputReferenceHarness {
    static final int NX = 4, NY = 3, NZ = 5, CLASSES = 6, APPROX = 4;
    static final String[] NAMES = {
        "Thalamus", "Cerebrum-GM", "thalamus", "Thalamus ", "Thalamus", "Background"
    };

    static float membership(int x, int y, int z, int k) {
        if ((x + 2*y + z + k) % 11 == 0) return 0.0f;
        return (x*71 + y*19 + z*7 + k*11 + 1) / 100.0f;
    }

    static int best(int x, int y, int z, int slot) {
        int voxel = (x * NY + y) * NZ + z;
        int[][] patterns = {
            {0, 4, 1, 4}, // both duplicate classes; last candidate slot wins
            {4, 4, 3, 2}, // duplicate class with two matching candidate slots
            {1, 2, 3, 1}, // no exact Thalamus class selected
            {0, 2, 3, 1}  // first exact Thalamus class only
        };
        return patterns[voxel % patterns.length][slot];
    }

    static float[][][][] run(boolean useApprox) {
        float[][][][] memberships = new float[NX][NY][NZ][CLASSES];
        int[][][][] best = new int[NX][NY][NZ][APPROX];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            for (int k=0;k<CLASSES;k++) memberships[x][y][z][k] = membership(x,y,z,k);
            for (int a=0;a<APPROX;a++) best[x][y][z][a] = best(x,y,z,a);
        }
        float[][][][] out = new float[1][NX][NY][NZ];
        // Excerpt of the selected method, preserving both loop orders.
        for (int k=0;k<CLASSES;k++) {
            if (NAMES[k].equals("Thalamus")) {
                for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
                    if (useApprox) {
                        for (int a=0;a<APPROX;a++)
                            if (best[x][y][z][a]==k)
                                out[0][x][y][z] = memberships[x][y][z][a];
                    } else {
                        out[0][x][y][z] = memberships[x][y][z][k];
                    }
                }
            }
        }
        return out;
    }

    static void floats(float[][][][] output) {
        System.out.print("[");
        boolean first = true;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            if (!first) System.out.print(",");
            System.out.print(Float.toString(output[0][x][y][z]));
            first = false;
        }
        System.out.print("]");
    }

    public static void main(String[] args) {
        System.out.print("{\"shape_xyz\":["+NX+","+NY+","+NZ+"],\"atlas_names\":[");
        for (int k=0;k<NAMES.length;k++) {
            if (k>0) System.out.print(",");
            System.out.print("\""+NAMES[k]+"\"");
        }
        System.out.print("],\"approx_slots\":"+APPROX+",\"java_reference\":{\"exact\":");
        floats(run(false));
        System.out.print(",\"approx\":");
        floats(run(true));
        System.out.println("}}");
    }
}
