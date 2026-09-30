import java.util.Arrays;

/**
 * Small reference harness copied from LesionToadSegmentation's selected
 * generateDuraRemovalOutputs method. Inputs are deterministic synthetic
 * float[][][][] memberships and three companion maps; output is flattened in
 * [channel,x,y,z] order. See tests/data/synthetic/lesion_dura_outputs/README.md.
 */
public final class LesionToadsDuraRemovalOutputsReferenceHarness {
    private static final String[] names = {
        "SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM",
        "CerebralGM", "Cerebellum-GM", "CerebellarGM", "CerebellarGM)",
        "WhiteMatter", "Background", "Other", "other"
    };

    public static void main(String[] args) {
        final int nx = 5, ny = 2, nz = 2, classes = names.length;
        float[][][][] mems = new float[nx][ny][nz][classes];
        int[][][] segmentation = new int[nx][ny][nz];
        float[][][] lesionMembership = new float[nx][ny][nz];
        float[][][] blackHoleMem = new float[nx][ny][nz];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            segmentation[x][y][z] = (x + 2 * y + 3 * z) % classes;
            lesionMembership[x][y][z] = (float)(0.013 * (1 + x + 3*y + 7*z));
            blackHoleMem[x][y][z] = (float)(0.017 * (2 + 2*x + y + 5*z));
            for (int k = 0; k < classes; k++) {
                mems[x][y][z][k] = (float)(0.001 * (1 + x + 5*y + 11*z + 13*k));
            }
        }

        float[][][][] Mems = new float[4][nx][ny][nz];
        // CSF: the last matching alias overwrites the earlier one.
        for (int k = 0; k < classes; k++) {
            if (names[k].equals("SulcalCSF") || names[k].equals("Sulcal-CSF")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                    Mems[0][x][y][z] = mems[x][y][z][k];
            }
        }
        // Lesion-specific initialization and source class-exclusion loop.
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
            Mems[1][x][y][z] = 0.0f;
            Mems[2][x][y][z] = (lesionMembership[x][y][z] + blackHoleMem[x][y][z]);
            Mems[3][x][y][z] = 0.0f;
        }
        for (int k=0;k<classes;k++) {
            if (!names[k].equals("SulcalCSF")
                && !names[k].equals("Sulcal-CSF")
                && !names[k].equals("CorticalGM")
                && !names[k].equals("Cerebrum-GM")
                && !names[k].equals("CerebralGM")
                && !names[k].equals("Cerebellum-GM")
                && !names[k].equals("CerebellarGM")
                && !names[k].equals("Background")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    Mems[2][x][y][z] += mems[x][y][z][k];
                    if (segmentation[x][y][z] == k) Mems[3][x][y][z] = 1.0f;
                }
            }
        }
        // GM probabilities and the literal CerebellarGM) hard-mask predicate.
        for (int k=0;k<classes;k++) {
            if (names[k].equals("CorticalGM")
                || names[k].equals("Cerebrum-GM")
                || names[k].equals("CerebralGM")
                || names[k].equals("Cerebellum-GM")
                || names[k].equals("CerebellarGM")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    Mems[1][x][y][z] += mems[x][y][z][k];
                    if (names[segmentation[x][y][z]].equals("Cerebellum-GM")
                        || names[segmentation[x][y][z]].equals("CerebellarGM)"))
                        Mems[3][x][y][z] = 1.0f;
                }
            }
        }
        System.out.println("{\"source_method\":\"LesionToadSegmentation.generateDuraRemovalOutputs\","
            +"\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\","
            +"\"shape\":["+nx+","+ny+","+nz+"],\"names\":"+quotedNames()+",\"channels\":[");
        for (int c = 0; c < 4; c++) {
            if (c > 0) System.out.println(",");
            System.out.print("[");
            boolean first = true;
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                if (!first) System.out.print(",");
                System.out.print(Float.toString(Mems[c][x][y][z]));
                first = false;
            }
            System.out.print("]");
        }
        System.out.println("]}");
    }

    private static String quotedNames() {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < names.length; i++) {
            if (i > 0) out.append(',');
            out.append('"').append(names[i]).append('"');
        }
        return out.append(']').toString();
    }
}
