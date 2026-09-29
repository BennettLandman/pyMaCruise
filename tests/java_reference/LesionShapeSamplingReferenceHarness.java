import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/** Synthetic reference for LesionToadAtlas.getRegularTransformedShape. */
public final class LesionShapeSamplingReferenceHarness {
    private static final int[] IMAGE_DIMS = {7, 6, 5};
    private static final int[] TEMPLATE_DIMS = {6, 5, 4};
    private static final float[] IMAGE_RES = {1.0f, 1.25f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9375f, 0.9375f, 0.9375f};
    private static final int[] LOW = {0, 0, 0};
    private static final int[] HIGH = {7, 6, 5};

    private static final class Case {
        final String name;
        final int x, y, z, clazz;
        final int[] low, high;
        final float[] transform;
        Case(String name, int x, int y, int z, int clazz, int[] low, int[] high, float[] transform) {
            this.name = name; this.x = x; this.y = y; this.z = z; this.clazz = clazz;
            this.low = low; this.high = high; this.transform = transform;
        }
    }

    private static float[][][] shape(int clazz) {
        float[][][] values = new float[TEMPLATE_DIMS[0]][TEMPLATE_DIMS[1]][TEMPLATE_DIMS[2]];
        float offset = clazz == 0 ? 0.12f : 0.23f;
        for (int x = 0; x < TEMPLATE_DIMS[0]; x++) {
            for (int y = 0; y < TEMPLATE_DIMS[1]; y++) {
                for (int z = 0; z < TEMPLATE_DIMS[2]; z++) {
                    values[x][y][z] = offset + 0.017f*x + 0.023f*y + 0.031f*z + 0.004f*x*y;
                }
            }
        }
        return values;
    }

    private static float[] sample(Case c, float[][][] shape, float[] mappedOutput) {
        float defaultValue = c.clazz == 0 ? 1.0f : 0.0f;
        if (c.x < c.low[0] || c.x >= c.high[0]
                || c.y < c.low[1] || c.y >= c.high[1]
                || c.z < c.low[2] || c.z >= c.high[2]) {
            return new float[] {defaultValue, Float.NaN, Float.NaN, Float.NaN};
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
        float value = ImageFunctions.linearInterpolation(shape, defaultValue,
                mapped[0], mapped[1], mapped[2], TEMPLATE_DIMS[0], TEMPLATE_DIMS[1], TEMPLATE_DIMS[2]);
        mappedOutput[0] = mapped[0]; mappedOutput[1] = mapped[1]; mappedOutput[2] = mapped[2];
        return new float[] {value, mapped[0], mapped[1], mapped[2]};
    }

    private static String vector(int[] values) {
        return String.format("[%d,%d,%d]", values[0], values[1], values[2]);
    }

    private static String vector(float[] values) {
        return String.format("[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g]", values[0], values[1],
                values[2], values[3], values[4], values[5]);
    }

    public static void main(String[] args) {
        Case[] cases = {
            new Case("identity_interior", 3, 2, 2, 1, LOW, HIGH, new float[] {0,0,0,0,0,0}),
            new Case("rigid_anisotropic_class0", 3, 2, 2, 0, LOW, HIGH,
                    new float[] {0.12f,-0.08f,0.04f,0.35f,-0.2f,0.1f}),
            new Case("input_bbox_default_class0", 6, 4, 3, 0, LOW, new int[] {6,6,5},
                    new float[] {0,0,0,0,0,0}),
            new Case("input_bbox_default_class1", 0, 4, 3, 1, new int[] {1,0,0}, HIGH,
                    new float[] {0,0,0,0,0,0}),
            new Case("template_low_default_class0", 2, 2, 2, 0, LOW, HIGH,
                    new float[] {0,0,0,-30,-30,-30}),
            new Case("template_low_default_class1", 2, 2, 2, 1, LOW, HIGH,
                    new float[] {0,0,0,-30,-30,-30}),
            new Case("template_high_default_class1", 2, 2, 2, 1, LOW, HIGH,
                    new float[] {0,0,0,30,30,30}),
            new Case("rigid_fractional_class1", 4, 3, 2, 1, LOW, HIGH,
                    new float[] {0.07f,-0.03f,0.09f,0.11f,0.17f,-0.13f})
        };
        System.out.print("{\"source\":\"VandyCRUISE LesionToadAtlas regular shape wrapper + ParametricTransform + ImageFunctions\",");
        System.out.print("\"image_dims\":[7,6,5],\"template_dims\":[6,5,4],");
        System.out.print("\"image_resolutions\":[1.0,1.25,1.5],\"template_resolutions\":[0.9375,0.9375,0.9375],");
        System.out.print("\"shape_formula\":\"offset+0.017*x+0.023*y+0.031*z+0.004*x*y; offset=0.12(class0)/0.23(other)\",");
        System.out.print("\"cases\":[");
        for (int n = 0; n < cases.length; n++) {
            Case c = cases[n];
            float[] mapped = new float[3];
            float[] output = sample(c, shape(c.clazz), mapped);
            if (n > 0) System.out.print(",");
            String mappedJson = Float.isNaN(output[1]) ? "null"
                    : String.format("[%.9g,%.9g,%.9g]", output[1], output[2], output[3]);
            System.out.printf("{\"name\":\"%s\",\"point\":[%d,%d,%d],\"class\":%d,"
                            + "\"bbox_low\":%s,\"bbox_high\":%s,\"transform\":%s,"
                            + "\"mapped\":%s,\"value\":%.9g}",
                    c.name, c.x, c.y, c.z, c.clazz, vector(c.low), vector(c.high),
                    vector(c.transform), mappedJson, output[0]);
        }
        System.out.println("]}");
    }
}
