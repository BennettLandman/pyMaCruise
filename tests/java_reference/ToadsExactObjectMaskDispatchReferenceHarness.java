import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/**
 * Fixed-input excerpt of the OBJECT/MASK dispatch in selected
 * MultipleToadSegmentation.computeExactMemberships.
 *
 * The test phantom is loaded from the committed RAS NIfTI fixture. Other
 * statistics are deliberately explicit constants, so this harness checks
 * class-term dispatch, in-place smoothing, relation scaling, normalization,
 * and objective order without claiming variance estimation or full TOADS.
 */
public final class ToadsExactObjectMaskDispatchReferenceHarness {
    private static final String SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX = 15, NY = 13, NZ = 11, CHANNELS = 2, CLASSES = 2;
    private static final float SMOOTH = 0.03f, FACTOR = 0.62f, COEFFICIENT = 0.08f;
    private static final float ZERO = 1e-10f, INF = 1e10f;
    private static final float[] MODALITY = {0.7f, 0.3f};
    private static final float[] IMAGE_VARIANCE = {0.75f, 0.25f};
    private static final float[][] VARIANCE = {{0.8f, 0.9f}, {1.1f, 0.6f}};
    private static final float[][] CENTROID = {{0.42f, 0.57f}, {0.61f, 0.28f}};
    private static final float[] SMOOTH_PRIOR = {0.6f, 0.85f};
    private static final float[] ATLAS_PRIOR = {0.45f, 0.55f};
    private static final float[][] SIMILARITY = {{0.0f, 0.4f}, {0.5f, 0.0f}};
    private static float[][][][] images;
    private static int[][][] ribbon;
    private static float[][][][] memberships;
    private static float[][][][] shapes;
    private static int[][][][] closest;

    private static float[][][][] readNifti(String path) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (GZIPInputStream input = new GZIPInputStream(new FileInputStream(path))) {
            byte[] chunk = new byte[8192];
            int n;
            while ((n = input.read(chunk)) >= 0) bytes.write(chunk, 0, n);
        }
        byte[] data = bytes.toByteArray();
        ByteBuffer header = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);
        if (header.getInt(0) != 348 || header.getShort(70) != 16 || header.getShort(72) != 32)
            throw new IOException("fixture must be little-endian float32 NIfTI-1");
        if (header.getShort(40) != 4 || header.getShort(42) != NX || header.getShort(44) != NY ||
                header.getShort(46) != NZ || header.getShort(48) != CHANNELS)
            throw new IOException("unexpected fixture dimensions");
        int offset = (int) header.getFloat(108);
        float[][][][] values = new float[NX][NY][NZ][CHANNELS];
        ByteBuffer voxels = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);
        for (int c = 0; c < CHANNELS; c++) for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
                int index = (((c * NZ + z) * NY + y) * NX + x);
                values[x][y][z][c] = voxels.getFloat(offset + index * 4);
            }
        return values;
    }

    private static int[][][] readRibbon(String path) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (GZIPInputStream input = new GZIPInputStream(new FileInputStream(path))) {
            byte[] chunk = new byte[8192];
            int n;
            while ((n = input.read(chunk)) >= 0) bytes.write(chunk, 0, n);
        }
        byte[] data = bytes.toByteArray();
        ByteBuffer header = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);
        if (header.getInt(0) != 348 || header.getShort(70) != 2 || header.getShort(72) != 8 ||
                header.getShort(40) != 3 || header.getShort(42) != NX ||
                header.getShort(44) != NY || header.getShort(46) != NZ)
            throw new IOException("unexpected folded-ribbon NIfTI-1 fixture");
        int offset = (int) header.getFloat(108);
        int[][][] values = new int[NX][NY][NZ];
        for (int z = 0; z < NZ; z++) for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
            int index = (z * NY + y) * NX + x;
            values[x][y][z] = data[offset + index] & 0xff;
        }
        return values;
    }

    private static float shape(int x, int y, int z, int k) {
        float primary = 0.25f + 0.5f * ribbon[x][y][z];
        return k == 0 ? primary : 1.0f - primary;
    }

    private static float initial(int x, int y, int z, int k) {
        int n = k == 0 ? (x + y + z) % 5 : (x + 2 * y + z) % 7;
        return k == 0 ? 0.2f + n * 0.03f : 0.3f + n * 0.02f;
    }

    private static void initialize() {
        memberships = new float[NX][NY][NZ][CLASSES];
        shapes = new float[NX][NY][NZ][CLASSES];
        closest = new int[NX][NY][NZ][4];
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            for (int k = 0; k < CLASSES; k++) {
                memberships[x][y][z][k] = initial(x, y, z, k);
                shapes[x][y][z][k] = shape(x, y, z, k);
            }
            java.util.Arrays.fill(closest[x][y][z], -1);
            closest[x][y][z][0] = (x + y + z) % CLASSES;
            if ((x + 2 * y + z) % 3 == 0) closest[x][y][z][2] = 1 - closest[x][y][z][0];
        }
    }

    private static int[] update() {
        float total = 0.0f;
        int[][] offsets = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for (int x = 1; x < NX - 1; x++) for (int y = 1; y < NY - 1; y++) for (int z = 1; z < NZ - 1; z++) {
            float denominator = 0.0f;
            float[] objective = new float[CLASSES];
            for (int k = 0; k < CLASSES; k++) {
                float num = 0.0f;
                for (int c = 0; c < CHANNELS; c++) {
                    float weight = (MODALITY[c] * VARIANCE[c][k]) * IMAGE_VARIANCE[c];
                    if (k == 1) {
                        float distance = images[x][y][z][c] == 0.0f ? 0.0f : 1.0f;
                        num += weight * distance;
                    } else {
                        float residual = images[x][y][z][c] - CENTROID[c][k];
                        float weightedResidual = weight * residual;
                        num += weightedResidual * residual;
                    }
                }
                float energy = num;
                float neighbor = 0.0f;
                for (int[] d : offsets) for (int m = 0; m < CLASSES; m++) if (m != k) {
                    float value = memberships[x+d[0]][y+d[1]][z+d[2]][m];
                    neighbor += value * value;
                }
                if (SMOOTH > 0.0f) {
                    num += (SMOOTH_PRIOR[k] * SMOOTH) * neighbor;
                    energy += ((0.5f * SMOOTH_PRIOR[k]) * SMOOTH) * neighbor;
                }
                float inversePrior = 0.0f;
                for (int m = 0; m < CLASSES; m++) if (m != k) {
                    float term = ATLAS_PRIOR[m] * SIMILARITY[k][m];
                    term *= shapes[x][y][z][m];
                    term *= shapes[x][y][z][m];
                    inversePrior += term;
                }
                num += COEFFICIENT * inversePrior;
                energy += (0.5f * COEFFICIENT) * inversePrior;
                boolean direct = closest[x][y][z][0] == k || closest[x][y][z][1] == k;
                boolean neighborClass = closest[x][y][z][2] == k || closest[x][y][z][3] == k;
                float factor;
                if (direct) factor = 1.0f;
                else {
                    float scaled = (2.0f * FACTOR) * shapes[x][y][z][k];
                    factor = neighborClass ? Math.max(0.25f, Math.min(scaled, 1.0f)) : Math.min(1.0f, scaled);
                }
                float reciprocal = num > ZERO ? 1.0f / num : INF;
                memberships[x][y][z][k] = factor * reciprocal;
                denominator += memberships[x][y][z][k];
                objective[k] = energy;
            }
            if (denominator > 0.0f) for (int k = 0; k < CLASSES; k++) {
                memberships[x][y][z][k] = memberships[x][y][z][k] / denominator;
                float contribution = objective[k] * memberships[x][y][z][k];
                contribution *= memberships[x][y][z][k];
                total += contribution;
            }
        }
        int[] output = new int[NX * NY * NZ * CLASSES + 1];
        int index = 0;
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++)
            for (int k = 0; k < CLASSES; k++) output[index++] = Float.floatToIntBits(memberships[x][y][z][k]);
        output[index] = Float.floatToIntBits(total);
        return output;
    }

    public static void main(String[] args) throws Exception {
        images = readNifti(args[0]);
        ribbon = readRibbon(args[1]);
        initialize();
        int[] bits = update();
        System.out.print("{\"source_sha256\":\"" + SOURCE_SHA + "\",\"shape\":[15,13,11,2],\"class_types\":[\"OBJECT\",\"MASK\"],\"updated_bits\":[");
        for (int i = 0; i < bits.length - 1; i++) { if (i > 0) System.out.print(','); System.out.print(bits[i]); }
        System.out.println("],\"objective_bits\":" + bits[bits.length - 1] + "}");
    }
}
