import gov.nih.mipav.model.structures.jama.JamaMatrix;
import edu.jhmi.rad.medic.utilities.PowerTable;

import java.io.DataInputStream;
import java.io.FileInputStream;

/**
 * Isolated selected-source reference for one exact-membership, degree-one
 * global Chebyshev field fit.  The loop and arithmetic order follow
 * VandyCRUISE InhomogeneityCorrection2.computeExact3DPolynomialField and
 * AlgorithmMultipleToads' selected short-constructor defaults.  It is not an
 * end-to-end TOADS reference and does not exercise degree two or three.
 */
public final class ToadsPostThinningExactFieldReferenceHarness {
    private static int fbits(float value) { return Float.floatToRawIntBits(value); }

    private static String intArray(int[] a) {
        StringBuilder s = new StringBuilder("[");
        for (int i = 0; i < a.length; i++) { if (i > 0) s.append(','); s.append(a[i]); }
        return s.append(']').toString();
    }

    private static String doubleArray(double[] a) {
        StringBuilder s = new StringBuilder("[");
        for (int i = 0; i < a.length; i++) { if (i > 0) s.append(','); s.append(Double.toString(a[i])); }
        return s.append(']').toString();
    }

    private static String doubleMatrix(double[][] a) {
        StringBuilder s = new StringBuilder("[");
        for (int i = 0; i < a.length; i++) { if (i > 0) s.append(','); s.append(doubleArray(a[i])); }
        return s.append(']').toString();
    }

    public static void main(String[] args) throws Exception {
        int nx = Integer.parseInt(args[1]), ny = Integer.parseInt(args[2]), nz = Integer.parseInt(args[3]);
        int classes = Integer.parseInt(args[4]), volume = nx * ny * nz;
        byte[] labels = new byte[volume];
        float[] image = new float[volume], memberships = new float[volume * classes], field = new float[volume];
        float[] centers = new float[classes];
        try (DataInputStream in = new DataInputStream(new FileInputStream(args[0]))) {
            for (int i = 0; i < volume; i++) labels[i] = in.readByte();
            for (int i = 0; i < volume; i++) image[i] = in.readFloat();
            for (int i = 0; i < volume * classes; i++) memberships[i] = in.readFloat();
            for (int i = 0; i < volume; i++) field[i] = in.readFloat();
            for (int i = 0; i < classes; i++) centers[i] = in.readFloat();
        }

        double[][] matrixValues = new double[4][4];
        double[][] rhsValues = new double[4][1];
        JamaMatrix normal = new JamaMatrix(4, 4), rhs = new JamaMatrix(4, 1);
        // VandyCRUISE passes float literals to this double-parameter API.
        PowerTable csf = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        PowerTable gm = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        PowerTable wm = new PowerTable(0.0f, 2.0f, 0.000001f, 1.0f);
        for (int x = 0; x < nx; x += 3) for (int y = 0; y < ny; y += 3) for (int z = 0; z < nz; z += 3) {
            int voxel = (x * ny + y) * nz + z;
            if (labels[voxel] == 0) continue;
            float[] basis = {1.0f, (float)x / (float)(nx - 1), (float)y / (float)(ny - 1), (float)z / (float)(nz - 1)};
            for (int k = 1; k < classes; k++) {
                float corrected = field[voxel] * image[voxel];
                float residual = corrected - centers[k];
                float distance = residual * residual;
                float weight = 1.0f;
                if (distance > 1.0e-30f) {
                    PowerTable table = k <= 2 ? csf : (k <= 7 ? gm : wm);
                    weight = (float)(table.lookup((double)distance, 1.0f) / (double)distance);
                }
                float membershipSquared = memberships[voxel * classes + k] * memberships[voxel * classes + k];
                float vectorTerm = membershipSquared * weight;
                vectorTerm = vectorTerm * image[voxel];
                vectorTerm = vectorTerm * centers[k];
                float matrixTerm = membershipSquared * weight;
                matrixTerm = matrixTerm * image[voxel];
                matrixTerm = matrixTerm * image[voxel];
                for (int row = 0; row < 4; row++) {
                    rhsValues[row][0] += (double)(vectorTerm * basis[row]);
                    float rowTerm = matrixTerm * basis[row];
                    for (int col = 0; col < 4; col++) matrixValues[row][col] += (double)(rowTerm * basis[col]);
                }
            }
        }
        float fullVolume = (float)volume;
        normal = new JamaMatrix(4, 4);
        rhs = new JamaMatrix(4, 1);
        for (int row = 0; row < 4; row++) {
            rhs.set(row, 0, rhsValues[row][0] / (double)fullVolume);
            for (int col = 0; col < 4; col++) normal.set(row, col, matrixValues[row][col] / (double)fullVolume);
        }
        JamaMatrix coefficient = normal.solve(rhs);

        float mean = 0.0f, count = 0.0f;
        boolean[] mask = new boolean[volume];
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            int voxel = (x * ny + y) * nz + z;
            if (labels[voxel] == 0) continue;
            mask[voxel] = true;
            float[] basis = {1.0f, (float)x / (float)(nx - 1), (float)y / (float)(ny - 1), (float)z / (float)(nz - 1)};
            float value = 0.0f;
            for (int n = 0; n < 4; n++) value = (float)((double)coefficient.get(n, 0) * (double)basis[n] + (double)value);
            field[voxel] = value;
            mean = mean + value;
            count = count + 1.0f;
        }
        mean = mean / count;
        float offset = 1.0f - mean;
        for (int voxel = 0; voxel < volume; voxel++) if (mask[voxel]) field[voxel] = field[voxel] + offset;

        int[] fieldBits = new int[volume];
        double[] coefficients = new double[4];
        for (int i = 0; i < 4; i++) coefficients[i] = coefficient.get(i, 0);
        for (int i = 0; i < volume; i++) fieldBits[i] = fbits(field[i]);
        double[][] normalCopy = normal.getArrayCopy();
        double[][] rhsCopy = rhs.getArrayCopy();
        for (int i = 0; i < 4; i++) matrixValues[i] = normalCopy[i].clone();
        for (int i = 0; i < 4; i++) rhsValues[i][0] = rhsCopy[i][0];
        System.out.println("{\"reference_kind\":\"formula_harness_not_source_class_execution\",\"inhomogeneity_upstream_source_sha256\":\"6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8\",\"wrapper_upstream_source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"compiled_jama_matrix_copy_sha256\":\"76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c\",\"compiled_jama_lu_copy_sha256\":\"9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6\",\"compiled_jamamat_wrapper_copy_sha256\":\"efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526\",\"compiled_power_table_copy_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"vandy_power_table_upstream_source_sha256\":\"5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678\",\"shape\":[" + nx + "," + ny + "," + nz + "],\"degree\":1,\"stride\":3,\"background\":0,\"normal\":" + doubleMatrix(matrixValues) + ",\"rhs\":" + doubleArray(new double[]{rhsValues[0][0],rhsValues[1][0],rhsValues[2][0],rhsValues[3][0]}) + ",\"coefficients\":" + doubleArray(coefficients) + ",\"field_bits\":" + intArray(fieldBits) + "}");
    }
}
