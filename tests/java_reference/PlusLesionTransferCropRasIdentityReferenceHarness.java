import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * RAS identity-grid reference for the selected transfer/crop wrapper slice.
 *
 * Source anchors: ToadsCruiseSpaceToOriginalSpace calls
 * ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR), and CropBasedOnGM calls
 * CubicVolumeCropper.crop(GM, 0, edge), then reuses the resulting parameters
 * for each target. On an identical RAS grid, the source reorientation matrix
 * is identity and nearest-neighbor transfer preserves each voxel. This
 * harness does not model changed origins, shapes, or spacings and is not a
 * compiled execution of the dependency-bearing JIST classes.
 *
 * Input: nx,ny,nz,edge, then eight float32 target volumes in x/y/z order;
 * target 0 is multi-atlas GM and defines the crop. Output: low/high bounds,
 * crop dimensions, then all eight cropped volumes in the same order.
 */
public final class PlusLesionTransferCropRasIdentityReferenceHarness {
    private static int index(int x, int y, int z, int ny, int nz) {
        return (x * ny + y) * nz + z;
    }

    public static void main(String[] args) throws Exception {
        try (DataInputStream in = new DataInputStream(new FileInputStream(args[0]));
             DataOutputStream out = new DataOutputStream(new FileOutputStream(args[1]))) {
            int nx = in.readInt(), ny = in.readInt(), nz = in.readInt(), edge = in.readInt();
            int voxels = nx * ny * nz;
            float[][] volumes = new float[8][voxels];
            for (int v = 0; v < volumes.length; v++)
                for (int i = 0; i < voxels; i++) volumes[v][i] = in.readFloat();

            int[] low = {nx, ny, nz}, high = {-1, -1, -1};
            for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
                if (volumes[0][index(x, y, z, ny, nz)] > 0.0f) {
                    low[0] = Math.min(low[0], x); high[0] = Math.max(high[0], x);
                    low[1] = Math.min(low[1], y); high[1] = Math.max(high[1], y);
                    low[2] = Math.min(low[2], z); high[2] = Math.max(high[2], z);
                }
            }
            for (int a = 0; a < 3; a++) { low[a] -= edge; high[a] += edge; }
            int cx = high[0] - low[0] + 1, cy = high[1] - low[1] + 1, cz = high[2] - low[2] + 1;
            out.writeInt(low[0]); out.writeInt(high[0]);
            out.writeInt(low[1]); out.writeInt(high[1]);
            out.writeInt(low[2]); out.writeInt(high[2]);
            out.writeInt(cx); out.writeInt(cy); out.writeInt(cz);
            for (float[] volume : volumes) {
                for (int x = 0; x < cx; x++) for (int y = 0; y < cy; y++) for (int z = 0; z < cz; z++) {
                    int sx = x + low[0], sy = y + low[1], sz = z + low[2];
                    float value = 0.0f;
                    if (sx >= 0 && sx < nx && sy >= 0 && sy < ny && sz >= 0 && sz < nz)
                        value = volume[index(sx, sy, sz, ny, nz)];
                    out.writeFloat(value);
                }
            }
        }
    }
}
