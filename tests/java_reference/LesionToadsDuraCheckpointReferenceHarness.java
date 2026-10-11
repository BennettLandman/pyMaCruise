import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Applies the selected LesionToadSegmentation.generateDuraRemovalOutputs
 * loops to arrays read from an exact persisted Python workflow checkpoint.
 * Pinned source: VandyCRUISE LesionToadSegmentation.java, method
 * generateDuraRemovalOutputs, SHA-256
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 *
 * Binary input: big-endian nx,ny,nz,classes; memberships in [x,y,z,k]
 * order; int32 segmentation [x,y,z]; lesion and black-hole maps [x,y,z].
 * Remaining arguments are the exact source atlas names in class order.
 * Binary output is four big-endian float32 channels in [c,x,y,z] order.
 */
public final class LesionToadsDuraCheckpointReferenceHarness {
    public static void main(String[] args) throws Exception {
        if (args.length < 4) throw new IllegalArgumentException("input output and atlas names required");
        final int nx, ny, nz, classes;
        final int voxels;
        float[][][][] mems;
        int[][][] segmentation;
        float[][][] lesionMembership, blackHoleMem;
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(args[0])))) {
            nx = in.readInt(); ny = in.readInt(); nz = in.readInt(); classes = in.readInt();
            if (nx <= 0 || ny <= 0 || nz <= 0 || classes != args.length - 2)
                throw new IllegalArgumentException("invalid grid or atlas-name count");
            voxels = Math.multiplyExact(Math.multiplyExact(nx, ny), nz);
            mems = new float[nx][ny][nz][classes];
            segmentation = new int[nx][ny][nz];
            lesionMembership = new float[nx][ny][nz];
            blackHoleMem = new float[nx][ny][nz];
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                for (int k=0;k<classes;k++) mems[x][y][z][k] = in.readFloat();
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) segmentation[x][y][z] = in.readInt();
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) lesionMembership[x][y][z] = in.readFloat();
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) blackHoleMem[x][y][z] = in.readFloat();
            if (in.read() != -1) throw new IllegalArgumentException("trailing bytes in input");
        }

        final String[] names = new String[classes];
        System.arraycopy(args, 2, names, 0, classes);
        float[][][][] out = new float[4][nx][ny][nz];
        // Source: LesionToadSegmentation.generateDuraRemovalOutputs().
        for (int k=0;k<classes;k++) {
            if (names[k].equals("SulcalCSF") || names[k].equals("Sulcal-CSF"))
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                    out[0][x][y][z] = mems[x][y][z][k];
        }
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
            out[1][x][y][z] = 0.0f;
            out[2][x][y][z] = lesionMembership[x][y][z] + blackHoleMem[x][y][z];
            out[3][x][y][z] = 0.0f;
        }
        for (int k=0;k<classes;k++) {
            if (!names[k].equals("SulcalCSF") && !names[k].equals("Sulcal-CSF")
                && !names[k].equals("CorticalGM") && !names[k].equals("Cerebrum-GM")
                && !names[k].equals("CerebralGM") && !names[k].equals("Cerebellum-GM")
                && !names[k].equals("CerebellarGM") && !names[k].equals("Background")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    out[2][x][y][z] += mems[x][y][z][k];
                    if (segmentation[x][y][z] == k) out[3][x][y][z] = 1.0f;
                }
            }
        }
        for (int k=0;k<classes;k++) {
            if (names[k].equals("CorticalGM") || names[k].equals("Cerebrum-GM")
                || names[k].equals("CerebralGM") || names[k].equals("Cerebellum-GM")
                || names[k].equals("CerebellarGM")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    out[1][x][y][z] += mems[x][y][z][k];
                    if (names[segmentation[x][y][z]].equals("Cerebellum-GM")
                        || names[segmentation[x][y][z]].equals("CerebellarGM)"))
                        out[3][x][y][z] = 1.0f;
                }
            }
        }
        try (DataOutputStream output = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(args[1])))) {
            for (int c=0;c<4;c++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                output.writeFloat(out[c][x][y][z]);
        }
    }
}
