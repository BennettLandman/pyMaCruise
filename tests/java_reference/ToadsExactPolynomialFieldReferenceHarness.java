import gov.nih.mipav.model.structures.jama.JamaMatrix;
import edu.jhmi.rad.medic.utilities.PowerTable;

import java.io.DataInputStream;
import java.io.FileInputStream;

/** Selected-source formula harness for exact-membership TOADS field degrees 1-3. */
public final class ToadsExactPolynomialFieldReferenceHarness {
    private static int bits(float value) { return Float.floatToRawIntBits(value); }
    private static String floats(float[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(bits(values[i])); }
        return out.append(']').toString();
    }
    private static String ints(int[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(values[i]); }
        return out.append(']').toString();
    }
    private static String doubles(double[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(Double.toString(values[i])); }
        return out.append(']').toString();
    }
    private static float[] basis(int x, int y, int z, int nx, int ny, int nz, int degree) {
        float rx = (float)x / (float)(nx - 1), ry = (float)y / (float)(ny - 1), rz = (float)z / (float)(nz - 1);
        float[] p = new float[degree == 1 ? 4 : (degree == 2 ? 10 : 20)];
        p[0] = 1.0f; p[1] = rx; p[2] = ry; p[3] = rz;
        if (degree >= 2) {
            p[4] = 2.0f * rx * rx - 1.0f; p[5] = 2.0f * ry * ry - 1.0f; p[6] = 2.0f * rz * rz - 1.0f;
            p[7] = p[1] * p[2]; p[8] = p[2] * p[3]; p[9] = p[3] * p[1];
        }
        if (degree >= 3) {
            p[10] = 4.0f * rx * rx * rx - 3.0f * rx; p[11] = 4.0f * ry * ry * ry - 3.0f * ry; p[12] = 4.0f * rz * rz * rz - 3.0f * rz;
            p[13] = p[4] * p[2]; p[14] = p[1] * p[5]; p[15] = p[5] * p[3]; p[16] = p[2] * p[6]; p[17] = p[6] * p[1]; p[18] = p[3] * p[4]; p[19] = p[1] * p[2] * p[3];
        }
        return p;
    }
    public static void main(String[] args) throws Exception {
        int nx = Integer.parseInt(args[1]), ny = Integer.parseInt(args[2]), nz = Integer.parseInt(args[3]);
        int classes = Integer.parseInt(args[4]), degree = Integer.parseInt(args[5]), volume = nx * ny * nz;
        int np = degree == 1 ? 4 : (degree == 2 ? 10 : 20);
        byte[] labels = new byte[volume]; float[] image = new float[volume], mem = new float[volume * classes];
        float[] field = new float[volume], centers = new float[classes];
        try (DataInputStream input = new DataInputStream(new FileInputStream(args[0]))) {
            for (int i = 0; i < volume; i++) labels[i] = input.readByte();
            for (int i = 0; i < volume; i++) image[i] = input.readFloat();
            for (int i = 0; i < volume * classes; i++) mem[i] = input.readFloat();
            for (int i = 0; i < volume; i++) field[i] = input.readFloat();
            for (int i = 0; i < classes; i++) centers[i] = input.readFloat();
        }
        PowerTable csf = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        PowerTable gm = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        PowerTable wm = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        JamaMatrix normal = new JamaMatrix(np, np, 0.0f), rhs = new JamaMatrix(np, 1, 0.0f);
        for (int x = 0; x < nx; x += 3) for (int y = 0; y < ny; y += 3) for (int z = 0; z < nz; z += 3) {
            int voxel = (x * ny + y) * nz + z;
            if (labels[voxel] == 0) continue;
            float[] p = basis(x, y, z, nx, ny, nz, degree);
            for (int k = 1; k < classes; k++) {
                float corrected = field[voxel] * image[voxel], residual = corrected - centers[k];
                float distance = residual * residual, weight = 1.0f;
                if (distance > 1.0e-30f) {
                    PowerTable table = k <= 2 ? csf : (k <= 7 ? gm : wm);
                    weight = (float)(table.lookup((double)distance, 1.0f) / (double)distance);
                }
                float membership = mem[voxel * classes + k];
                float membershipSquared = membership * membership;
                float vectorTerm = membershipSquared * weight; vectorTerm = vectorTerm * image[voxel]; vectorTerm = vectorTerm * centers[k];
                float matrixTerm = membershipSquared * weight; matrixTerm = matrixTerm * image[voxel]; matrixTerm = matrixTerm * image[voxel];
                for (int row = 0; row < np; row++) {
                    rhs.set(row, 0, rhs.get(row, 0) + (double)(vectorTerm * p[row]));
                    float rowTerm = matrixTerm * p[row];
                    for (int col = 0; col < np; col++) normal.set(row, col, normal.get(row, col) + (double)(rowTerm * p[col]));
                }
            }
        }
        float fullVolume = (float)volume;
        for (int row = 0; row < np; row++) {
            rhs.set(row, 0, rhs.get(row, 0) / (double)fullVolume);
            for (int col = 0; col < np; col++) normal.set(row, col, normal.get(row, col) / (double)fullVolume);
        }
        JamaMatrix coefficient = normal.solve(rhs);
        float mean = 0.0f, count = 0.0f;
        boolean[] mask = new boolean[volume];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            int voxel = (x * ny + y) * nz + z; if (labels[voxel] == 0) continue; mask[voxel] = true;
            float[] p = basis(x, y, z, nx, ny, nz, degree); float value = 0.0f;
            for (int n = 0; n < np; n++) value = (float)((double)coefficient.get(n, 0) * (double)p[n] + (double)value);
            field[voxel] = value; mean = mean + value; count = count + 1.0f;
        }
        mean = mean / count; float offset = 1.0f - mean;
        for (int i = 0; i < volume; i++) if (mask[i]) field[i] = field[i] + offset;
        int[] fieldBits = new int[volume]; for (int i = 0; i < volume; i++) fieldBits[i] = bits(field[i]);
        double[] coefficients = new double[np]; for (int i = 0; i < np; i++) coefficients[i] = coefficient.get(i, 0);
        double[] normalFlat = new double[np * np], rhsFlat = new double[np];
        for (int i = 0; i < np; i++) { rhsFlat[i] = rhs.get(i, 0); for (int j = 0; j < np; j++) normalFlat[i * np + j] = normal.get(i, j); }
        System.out.println("{\"degree\":" + degree + ",\"stride\":3,\"normal\":" + doubles(normalFlat) + ",\"rhs\":" + doubles(rhsFlat) + ",\"coefficients\":" + doubles(coefficients) + ",\"field_bits\":" + ints(fieldBits) + "}");
    }
}
