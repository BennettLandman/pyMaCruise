import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Calls the selected VandyCRUISE transform/interpolation APIs on synthetic data. */
public final class ToadsShapeSamplingReferenceHarness {
    private static final int[] IMAGE_DIMS = {7, 6, 5};
    private static final int[] TEMPLATE_DIMS = {6, 5, 4};
    private static final float[] IMAGE_RES = {1.0f, 1.25f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9375f, 0.9375f, 0.9375f};

    private static final class Case {
        final String name;
        final int x, y, z, clazz;
        final int[] low, high;
        final float[] transform;
        final float slope;
        Case(String name, int x, int y, int z, int clazz,
                int[] low, int[] high, float[] transform, float slope) {
            this.name = name; this.x = x; this.y = y; this.z = z; this.clazz = clazz;
            this.low = low; this.high = high; this.transform = transform; this.slope = slope;
        }
    }

    private static float[][][] shape() {
        float[][][] values = new float[TEMPLATE_DIMS[0]][TEMPLATE_DIMS[1]][TEMPLATE_DIMS[2]];
        for (int x = 0; x < TEMPLATE_DIMS[0]; x++) {
            for (int y = 0; y < TEMPLATE_DIMS[1]; y++) {
                for (int z = 0; z < TEMPLATE_DIMS[2]; z++) {
                    values[x][y][z] = 0.2f + 0.02f*x + 0.025f*y + 0.03f*z + 0.005f*x*y;
                }
            }
        }
        return values;
    }

    private static float[] sample(Case c, float[][][] shape, float[] xyz) {
        if (c.x < c.low[0] || c.x >= c.high[0]
                || c.y < c.low[1] || c.y >= c.high[1]
                || c.z < c.low[2] || c.z >= c.high[2]) {
            return new float[] { c.clazz == 0 ? 1.0f : 0.0f, Float.NaN, Float.NaN, Float.NaN };
        }
        float x0i = IMAGE_DIMS[0] / 2.0f, y0i = IMAGE_DIMS[1] / 2.0f, z0i = IMAGE_DIMS[2] / 2.0f;
        float x0t = TEMPLATE_DIMS[0] / 2.0f, y0t = TEMPLATE_DIMS[1] / 2.0f, z0t = TEMPLATE_DIMS[2] / 2.0f;
        ParametricTransform model = new ParametricTransform("rigid",
                x0i, y0i, z0i, IMAGE_RES[0], IMAGE_RES[1], IMAGE_RES[2],
                IMAGE_DIMS[0], IMAGE_DIMS[1], IMAGE_DIMS[2],
                x0t, y0t, z0t, TEMPLATE_RES[0], TEMPLATE_RES[1], TEMPLATE_RES[2],
                TEMPLATE_DIMS[0], TEMPLATE_DIMS[1], TEMPLATE_DIMS[2]);
        float[][] rotation = ParametricTransform.computeRotation(c.transform);
        float[] mapped = model.imageToTemplate(c.x, c.y, c.z, c.transform, rotation, 1.0f);
        float interpolated = ImageFunctions.linearClosestInterpolation(
                shape, mapped[0], mapped[1], mapped[2],
                TEMPLATE_DIMS[0], TEMPLATE_DIMS[1], TEMPLATE_DIMS[2]);
        float result = Numerics.bounded(c.slope * (interpolated - 0.5f) + 0.5f, 0.0f, 1.0f);
        xyz[0] = mapped[0]; xyz[1] = mapped[1]; xyz[2] = mapped[2];
        return new float[] {result, mapped[0], mapped[1], mapped[2]};
    }

    private static String vector(int[] values) {
        return String.format("[%d,%d,%d]", values[0], values[1], values[2]);
    }

    private static String vector(float[] values) {
        return String.format("[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g]", values[0], values[1],
                values[2], values[3], values[4], values[5]);
    }

    public static void main(String[] args) {
        int[] low = {0, 0, 0}, high = {7, 6, 5};
        Case[] cases = {
            new Case("identity", 3, 2, 2, 0, low, high,
                    new float[] {0,0,0,0,0,0}, 2.0f),
            new Case("rigid_anisotropic", 4, 3, 2, 1, low, high,
                    new float[] {0.12f,-0.08f,0.04f,0.35f,-0.2f,0.1f}, 2.0f),
            new Case("closest_low_edge", 2, 2, 2, 1, low, high,
                    new float[] {0,0,0,-30,-30,-30}, 2.0f),
            new Case("closest_high_edge", 2, 2, 2, 1, low, high,
                    new float[] {0,0,0,30,30,30}, 2.0f),
            new Case("slope_low_clip", 0, 0, 0, 1, low, high,
                    new float[] {0,0,0,0,0,0}, 8.0f),
            new Case("slope_high_clip", 6, 5, 4, 1, low, high,
                    new float[] {0,0,0,0,0,0}, 8.0f),
            new Case("class0_bbox_default", 6, 4, 3, 0, low, new int[] {6,6,5},
                    new float[] {0,0,0,0,0,0}, 2.0f),
            new Case("class1_bbox_default", 0, 4, 3, 1, new int[] {1,0,0}, high,
                    new float[] {0,0,0,0,0,0}, 2.0f)
        };
        float[][][] shape = shape();
        System.out.print("{\"source\":\"VandyCRUISE ParametricTransform + ImageFunctions + shape postprocess\",");
        System.out.print("\"image_dims\":[7,6,5],\"template_dims\":[6,5,4],");
        System.out.print("\"image_resolutions\":[1.0,1.25,1.5],\"template_resolutions\":[0.9375,0.9375,0.9375],");
        System.out.print("\"shape\":\"0.2+0.02*x+0.025*y+0.03*z+0.005*x*y\",\"cases\":[");
        for (int n = 0; n < cases.length; n++) {
            Case c = cases[n];
            float[] mapped = new float[3];
            float[] output = sample(c, shape, mapped);
            if (n > 0) System.out.print(",");
            String mappedJson = Float.isNaN(output[1])
                    ? "null"
                    : String.format("[%.9g,%.9g,%.9g]", output[1], output[2], output[3]);
            System.out.printf("{\"name\":\"%s\",\"point\":[%d,%d,%d],\"class\":%d,"
                    + "\"bbox_low\":%s,\"bbox_high\":%s,\"transform\":%s,"
                    + "\"slope\":%.9g,\"mapped\":%s,\"value\":%.9g}",
                    c.name, c.x, c.y, c.z, c.clazz, vector(c.low), vector(c.high),
                    vector(c.transform), c.slope, mappedJson, output[0]);
        }
        System.out.println("]}");
    }
}
