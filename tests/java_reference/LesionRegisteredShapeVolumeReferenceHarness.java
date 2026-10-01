import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/**
 * Hand-coded Java excerpt for LesionToadAtlas.computeTransformedShapeBoundingBox,
 * getFastTransformedShape, and fastImageToShapeCoordinates. It uses the actual
 * pinned ParametricTransform and ImageFunctions Java sources for matrix setup
 * and interpolation; box values are supplied from the already Java-compared
 * registerShapes handoff. This does not execute the whole LesionToadAtlas class.
 *
 * Pinned source SHA-256 values:
 * LesionToadAtlas.java 718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e
 * ParametricTransform.java 22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac
 * ImageFunctions.java be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892
 */
public final class LesionRegisteredShapeVolumeReferenceHarness {
    private LesionRegisteredShapeVolumeReferenceHarness() {}

    public static void main(String[] args) throws Exception {
        try (DataInputStream in = new DataInputStream(new FileInputStream(args[0]));
             DataOutputStream out = new DataOutputStream(new FileOutputStream(args[1]))) {
            int nx = in.readInt(), ny = in.readInt(), nz = in.readInt(), classes = in.readInt();
            float rix = in.readFloat(), riy = in.readFloat(), riz = in.readFloat();
            float rtx = in.readFloat(), rty = in.readFloat(), rtz = in.readFloat();
            float[] transform = new float[6];
            for (int p = 0; p < transform.length; p++) transform[p] = in.readFloat();
            int[][] lower = new int[classes][3], upper = new int[classes][3];
            for (int k = 0; k < classes; k++) {
                for (int a = 0; a < 3; a++) lower[k][a] = in.readInt();
                for (int a = 0; a < 3; a++) upper[k][a] = in.readInt();
            }
            float[][][][] shapes = new float[classes][nx][ny][nz];
            for (int k = 0; k < classes; k++)
                for (int x = 0; x < nx; x++)
                    for (int y = 0; y < ny; y++)
                        for (int z = 0; z < nz; z++) shapes[k][x][y][z] = in.readFloat();

            ParametricTransform model = new ParametricTransform(
                "rigid", nx / 2.0f, ny / 2.0f, nz / 2.0f, rix, riy, riz, nx, ny, nz,
                shapes[0].length / 2.0f, shapes[0][0].length / 2.0f,
                shapes[0][0][0].length / 2.0f, rtx, rty, rtz,
                shapes[0].length, shapes[0][0].length, shapes[0][0][0].length
            );
            float[][] matrix = new float[3][4];
            float[][] rotation = ParametricTransform.computeRotation(transform);
            model.precomputeImageToTemplateMatrix(matrix, transform, rotation, 1.0f);

            for (int x = 0; x < nx; x++) {
                for (int y = 0; y < ny; y++) {
                    for (int z = 0; z < nz; z++) {
                        for (int k = 0; k < classes; k++) {
                            float value;
                            if (x < lower[k][0] || x >= upper[k][0]
                                || y < lower[k][1] || y >= upper[k][1]
                                || z < lower[k][2] || z >= upper[k][2]) {
                                value = k == 0 ? 1.0f : 0.0f;
                            } else {
                                // Exact fastImageToShapeCoordinates expression order.
                                float xp = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3];
                                float yp = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3];
                                float zp = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3];
                                value = ImageFunctions.linearInterpolation(
                                    shapes[k], k == 0 ? 1.0f : 0.0f, xp, yp, zp,
                                    shapes[k].length, shapes[k][0].length, shapes[k][0][0].length
                                );
                            }
                            out.writeFloat(value);
                        }
                    }
                }
            }
        }
    }
}
