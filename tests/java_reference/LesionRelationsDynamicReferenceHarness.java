import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/** Selected LesionToadSegmentation relation refresh on caller-supplied labels. */
public final class LesionRelationsDynamicReferenceHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("input.bin output.bin required");
        int nx, ny, nz, classes, limit;
        byte[][][] segmentation;
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(args[0])))) {
            nx = in.readInt(); ny = in.readInt(); nz = in.readInt(); classes = in.readInt(); limit = in.readInt();
            segmentation = new byte[nx][ny][nz];
            for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
                segmentation[x][y][z] = in.readByte();
            if (in.read() != -1) throw new IllegalArgumentException("unexpected trailing input");
        }
        LesionRelationsReferenceHarness.State state =
            LesionRelationsReferenceHarness.computeRelations(segmentation, classes, limit);
        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(args[1])))) {
            for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++)
                out.writeShort(state.map[x][y][z]);
            out.writeInt(state.list.length);
            for (boolean[] row : state.list) for (boolean related : row) out.writeByte(related ? 1 : 0);
            out.writeInt(state.iterations);
        }
    }
}
