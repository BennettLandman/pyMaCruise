import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Selected-source rigid mapping at landmarks on an RAS axis-aligned grid. */
public final class ToadsOrientationMappingHarness {
    private static String point(ParametricTransform model, RotationMatrix rotation,
                                int x, int y, int z, float[] transform) {
        float[] value = model.imageToTemplate(x, y, z, transform,
                rotation.getMatrix(), 1.0f);
        return "[" + value[0] + "," + value[1] + "," + value[2] + "]";
    }

    public static void main(String[] args) {
        int nx = 7, ny = 6, nz = 5;
        int ntx = 11, nty = 10, ntz = 9;
        // MultipleToadPreprocess AXIAL signs for RAS axis codes L2R/P2A/I2S:
        // signed resolutions are (-sx,-sy,+sz).
        ParametricTransform model = new ParametricTransform("rigid",
                nx / 2.0f, ny / 2.0f, nz / 2.0f,
                -1.25f, -1.75f, 2.25f, nx, ny, nz,
                ntx / 2.0f, nty / 2.0f, ntz / 2.0f,
                0.8f, 0.9f, 1.1f, ntx, nty, ntz);
        float[] transform = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        RotationMatrix rotation = model.computeRotationMatrix(transform);
        System.out.println("{\"signed_resolutions\":[-1.25,-1.75,2.25],"
                + "\"points\":["
                + point(model, rotation, 0, 0, 0, transform) + ","
                + point(model, rotation, 3, 2, 1, transform) + ","
                + point(model, rotation, 6, 5, 4, transform) + "]}");
    }
}
