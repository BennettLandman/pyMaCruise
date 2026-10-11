import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/**
 * Selected getFastTransformedShape / fastImageToShapeCoordinates excerpt for
 * P41z, where the registered image grid and template atlas grid differ.
 *
 * Input: image xyz dimensions, template xyz dimensions, class count, image and
 * template resolutions, six rigid parameters, per-class lower/upper bounds,
 * then each class's template-grid float map. Output is image-grid class maps
 * in x/y/z/class order. Matrix setup and interpolation use the pinned Java
 * ParametricTransform and ImageFunctions sources.
 */
public final class P41zRegisteredShapeVolumeReferenceHarness {
    private P41zRegisteredShapeVolumeReferenceHarness() {}

    public static void main(String[] args) throws Exception {
        try (DataInputStream in = new DataInputStream(new FileInputStream(args[0]));
             DataOutputStream out = new DataOutputStream(new FileOutputStream(args[1]))) {
            int ix = in.readInt(), iy = in.readInt(), iz = in.readInt();
            int tx = in.readInt(), ty = in.readInt(), tz = in.readInt();
            int classes = in.readInt();
            float rix = in.readFloat(), riy = in.readFloat(), riz = in.readFloat();
            float rtx = in.readFloat(), rty = in.readFloat(), rtz = in.readFloat();
            float[] transform = new float[6];
            for (int p = 0; p < transform.length; p++) transform[p] = in.readFloat();
            int[][] lower = new int[classes][3], upper = new int[classes][3];
            for (int k = 0; k < classes; k++) {
                for (int a = 0; a < 3; a++) lower[k][a] = in.readInt();
                for (int a = 0; a < 3; a++) upper[k][a] = in.readInt();
            }
            float[][][][] shapes = new float[classes][tx][ty][tz];
            for (int k = 0; k < classes; k++)
                for (int x = 0; x < tx; x++)
                    for (int y = 0; y < ty; y++)
                        for (int z = 0; z < tz; z++) shapes[k][x][y][z] = in.readFloat();
            if (in.read() != -1) throw new IllegalArgumentException("trailing input");

            ParametricTransform model = new ParametricTransform(
                "rigid", ix / 2.0f, iy / 2.0f, iz / 2.0f, rix, riy, riz, ix, iy, iz,
                tx / 2.0f, ty / 2.0f, tz / 2.0f, rtx, rty, rtz, tx, ty, tz
            );
            float[][] matrix = new float[3][4];
            float[][] rotation = ParametricTransform.computeRotation(transform);
            model.precomputeImageToTemplateMatrix(matrix, transform, rotation, 1.0f);

            for (int x = 0; x < ix; x++) {
                for (int y = 0; y < iy; y++) {
                    for (int z = 0; z < iz; z++) {
                        for (int k = 0; k < classes; k++) {
                            float value;
                            if (x < lower[k][0] || x >= upper[k][0]
                                || y < lower[k][1] || y >= upper[k][1]
                                || z < lower[k][2] || z >= upper[k][2]) {
                                value = k == 0 ? 1.0f : 0.0f;
                            } else {
                                // Preserve source float expression order.
                                float xp = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3];
                                float yp = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3];
                                float zp = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3];
                                value = ImageFunctions.linearInterpolation(
                                    shapes[k], k == 0 ? 1.0f : 0.0f, xp, yp, zp, tx, ty, tz
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
