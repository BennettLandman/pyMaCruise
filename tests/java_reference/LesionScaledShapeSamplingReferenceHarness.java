import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/** Selected Java scale-aware atlas sampling formulas for the lesion pyramid. */
public final class LesionScaledShapeSamplingReferenceHarness {
    private static final int[] IMAGE_DIMS = {21, 19, 17};
    private static final int[] TEMPLATE_DIMS = {24, 24, 24};
    private static final float[] IMAGE_RES = {1.1f, 1.3f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9f, 1.0f, 1.2f};
    private static final float[] TRANSFORM = {0.05f, -0.03f, 0.04f, 0.1f, 0.2f, -0.1f};

    private static String ints(int[] a) {
        return String.format("[%d,%d,%d]", a[0], a[1], a[2]);
    }

    private static String floats(float[] a) {
        return String.format("[%.9g,%.9g,%.9g]", a[0], a[1], a[2]);
    }

    private static String transform(float[] values) {
        return String.format("[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g]", values[0],
                values[1], values[2], values[3], values[4], values[5]);
    }

    private static float[][][] shape(int nx, int ny, int nz) {
        float[][][] image = new float[nx][ny][nz];
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++)
            image[x][y][z] = 0.1f + 0.01f*x + 0.02f*y + 0.03f*z + 0.001f*x*y;
        return image;
    }

    private static void appendCase(StringBuilder out, int scale, String name,
            int[] point, float[] trans) {
        int[] imageCoarse = {IMAGE_DIMS[0]/scale, IMAGE_DIMS[1]/scale, IMAGE_DIMS[2]/scale};
        int[] templateCoarse = {TEMPLATE_DIMS[0]/scale, TEMPLATE_DIMS[1]/scale,
                TEMPLATE_DIMS[2]/scale};
        float[] imageCenter = {IMAGE_DIMS[0]/2.0f, IMAGE_DIMS[1]/2.0f, IMAGE_DIMS[2]/2.0f};
        float[] templateCenter = {8.0f, 8.0f, 8.0f};
        ParametricTransform model = new ParametricTransform("rigid",
                imageCenter[0], imageCenter[1], imageCenter[2],
                IMAGE_RES[0], IMAGE_RES[1], IMAGE_RES[2],
                IMAGE_DIMS[0], IMAGE_DIMS[1], IMAGE_DIMS[2],
                templateCenter[0], templateCenter[1], templateCenter[2],
                TEMPLATE_RES[0], TEMPLATE_RES[1], TEMPLATE_RES[2],
                TEMPLATE_DIMS[0], TEMPLATE_DIMS[1], TEMPLATE_DIMS[2]);
        float[][] rotation = ParametricTransform.computeRotation(trans);
        float[] mapped = model.imageToTemplate(point[0], point[1], point[2], trans,
                rotation, (float)scale);
        float[][][] prior = shape(templateCoarse[0], templateCoarse[1], templateCoarse[2]);
        float value = ImageFunctions.linearInterpolation(prior, 0.0f,
                mapped[0], mapped[1], mapped[2], templateCoarse[0],
                templateCoarse[1], templateCoarse[2]);
        float[] derivative = {
            ImageFunctions.linearInterpolationXderivative(prior, mapped[0], mapped[1], mapped[2],
                    templateCoarse[0], templateCoarse[1], templateCoarse[2]),
            ImageFunctions.linearInterpolationYderivative(prior, mapped[0], mapped[1], mapped[2],
                    templateCoarse[0], templateCoarse[1], templateCoarse[2]),
            ImageFunctions.linearInterpolationZderivative(prior, mapped[0], mapped[1], mapped[2],
                    templateCoarse[0], templateCoarse[1], templateCoarse[2])};
        if (out.charAt(out.length()-1) != '[') out.append(',');
        out.append(String.format(
                "{\"scale\":%d,\"name\":\"%s\",\"point\":%s,\"image_dims\":%s,"
                + "\"template_dims\":%s,\"image_center\":%s,\"template_center\":%s,"
                + "\"image_resolutions\":[1.1,1.3,1.5],\"template_resolutions\":[0.9,1.0,1.2],"
                + "\"transform\":%s,\"mapped\":%s,\"value\":%.9g,\"derivative\":%s}",
                scale, name, ints(point), ints(imageCoarse), ints(templateCoarse),
                floats(imageCenter), floats(templateCenter), transform(trans), floats(mapped),
                value, floats(derivative)));
    }

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"cases\":[");
        for (int scale : new int[] {2, 4, 8}) {
            int[] coarse = {IMAGE_DIMS[0]/scale, IMAGE_DIMS[1]/scale, IMAGE_DIMS[2]/scale};
            int[] center = {coarse[0]/2, coarse[1]/2, coarse[2]/2};
            appendCase(out, scale, "fractional_interior", center, TRANSFORM);
            float[] low = TRANSFORM.clone(); low[3] = -100.0f;
            appendCase(out, scale, "low_constant_zero", center, low);
            float[] high = TRANSFORM.clone(); high[3] = 100.0f;
            appendCase(out, scale, "high_constant_zero", center, high);
        }
        System.out.println(out.append("]}"));
    }
}
