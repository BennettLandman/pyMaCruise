import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.lang.reflect.Field;
import edu.jhmi.rad.medic.methods.MultipleToadDeformableAtlas;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.libraries.ImageFunctions;

/** Executes the supplied VandyCRUISE pyramid registration entry point headlessly. */
public final class ToadsPyramidActualRegistrationHarness {
    private static float[] read(Path path, int count) throws Exception {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length != count * 4) throw new IllegalArgumentException("unexpected float payload: " + path);
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[count];
        for (int i = 0; i < count; i++) values[i] = buffer.getFloat();
        return values;
    }
    private static void set(Object target, String name, Object value) throws Exception {
        Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        field.set(target, value);
    }
    private static Object get(Object target, String name) throws Exception {
        Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.get(target);
    }
    private static String bits(float value) { return Integer.toString(Float.floatToIntBits(value)); }
    private static String vector(float[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) { if (i > 0) out.append(','); out.append(bits(values[i])); }
        return out.append(']').toString();
    }
    public static void main(String[] args) throws Exception {
        if (args.length != 23) throw new IllegalArgumentException(
            "usage: mem.bin shapes.bin nix niy niz ntx nty ntz classes rix riy riz rtx rty rtz levels iterations tx ty tz rx ry rz");
        int p = 0;
        Path memPath = Path.of(args[p++]), shapePath = Path.of(args[p++]);
        int nix = Integer.parseInt(args[p++]), niy = Integer.parseInt(args[p++]), niz = Integer.parseInt(args[p++]);
        int ntx = Integer.parseInt(args[p++]), nty = Integer.parseInt(args[p++]), ntz = Integer.parseInt(args[p++]);
        int classes = Integer.parseInt(args[p++]);
        float rix = Float.parseFloat(args[p++]), riy = Float.parseFloat(args[p++]), riz = Float.parseFloat(args[p++]);
        float rtx = Float.parseFloat(args[p++]), rty = Float.parseFloat(args[p++]), rtz = Float.parseFloat(args[p++]);
        int levels = Integer.parseInt(args[p++]), iterations = Integer.parseInt(args[p++]);
        float[] initial = new float[] { Float.parseFloat(args[p++]), Float.parseFloat(args[p++]), Float.parseFloat(args[p++]),
                                        Float.parseFloat(args[p++]), Float.parseFloat(args[p++]), Float.parseFloat(args[p++]) };
        float[] flatMem = read(memPath, nix * niy * niz * classes);
        float[] flatShapes = read(shapePath, ntx * nty * ntz * classes);
        float[][][][] mem = new float[nix][niy][niz][classes];
        float[][][][] shapes = new float[classes][ntx][nty][ntz];
        int q = 0;
        for (int x = 0; x < nix; x++) for (int y = 0; y < niy; y++) for (int z = 0; z < niz; z++) for (int k = 0; k < classes; k++)
            mem[x][y][z][k] = flatMem[q++];
        q = 0;
        for (int x = 0; x < ntx; x++) for (int y = 0; y < nty; y++) for (int z = 0; z < ntz; z++) for (int k = 0; k < classes; k++)
            shapes[k][x][y][z] = flatShapes[q++];

        MultipleToadDeformableAtlas atlas = new MultipleToadDeformableAtlas(classes);
        atlas.setImageInfo(nix, niy, niz, rix, riy, riz, 0, 0, 0, 0);
        set(atlas, "ntx", ntx); set(atlas, "nty", nty); set(atlas, "ntz", ntz);
        set(atlas, "rtx", rtx); set(atlas, "rty", rty); set(atlas, "rtz", rtz);
        set(atlas, "x0t", ntx / 2.0f); set(atlas, "y0t", nty / 2.0f); set(atlas, "z0t", ntz / 2.0f);
        set(atlas, "rsx", new float[classes]); set(atlas, "rsy", new float[classes]); set(atlas, "rsz", new float[classes]);
        atlas.setShapes(shapes);
        for (int k = 0; k < classes; k++) atlas.setShapeInfo(k, ntx, nty, ntz, rtx, rty, rtz);
        boolean[] registered = new boolean[classes];
        for (int k = 1; k < classes; k++) registered[k] = true;
        set(atlas, "registeredShape", registered);
        set(atlas, "transformModel", new ParametricTransform("rigid", nix / 2.0f, niy / 2.0f, niz / 2.0f,
            rix, riy, riz, nix, niy, niz, ntx / 2.0f, nty / 2.0f, ntz / 2.0f, rtx, rty, rtz, ntx, nty, ntz));
        set(atlas, "transformMode", 0);
        atlas.setTransform(initial);
        byte[][][][] best = new byte[nix][niy][niz][4];
        byte[][][] segmentation = new byte[nix][niy][niz];
        atlas.updateShapeRegistration(mem, best, segmentation, new String[0], "rigid", iterations, levels, 0.0f, 1.0f);
        atlas.setTransform(initial);
        atlas.registerShapesPyramid();
        float[][][][] coarseMem = ImageFunctions.subsample(mem, nix, niy, niz, classes, 4);
        float[][][][] coarseShapes = new float[classes][][][];
        int[] coarseX = new int[classes], coarseY = new int[classes], coarseZ = new int[classes];
        for (int k = 0; k < classes; k++) {
            coarseShapes[k] = ImageFunctions.subsample(shapes[k], ntx, nty, ntz, 4);
            coarseX[k] = ntx / 4; coarseY[k] = nty / 4; coarseZ[k] = ntz / 4;
        }
        java.lang.reflect.Method energyMethod = atlas.getClass().getDeclaredMethod(
            "computeExactScaledRegistrationEnergy", float[].class, float[][][][].class,
            float[][][][].class, int.class, int.class, int.class, int[].class, int[].class,
            int[].class, float.class);
        energyMethod.setAccessible(true);
        float[] nanTrial = new float[] { Float.NaN, Float.NaN, Float.NaN, Float.NaN, Float.NaN, Float.NaN };
        float coarseInitialEnergy = (Float)energyMethod.invoke(atlas, initial, coarseMem, coarseShapes,
            nix / 4, niy / 4, niz / 4, coarseX, coarseY, coarseZ, 4.0f);
        float coarseNanTrialEnergy = (Float)energyMethod.invoke(atlas, nanTrial, coarseMem, coarseShapes,
            nix / 4, niy / 4, niz / 4, coarseX, coarseY, coarseZ, 4.0f);
        float[][][][] fineMem = ImageFunctions.subsample(mem, nix, niy, niz, classes, 2);
        float[][][][] fineShapes = new float[classes][][][];
        int[] fineX = new int[classes], fineY = new int[classes], fineZ = new int[classes];
        for (int k = 0; k < classes; k++) {
            fineShapes[k] = ImageFunctions.subsample(shapes[k], ntx, nty, ntz, 2);
            fineX[k] = ntx / 2; fineY[k] = nty / 2; fineZ[k] = ntz / 2;
        }
        float finalEnergy = (Float)energyMethod.invoke(atlas, atlas.getTransform(), fineMem, fineShapes,
            nix / 2, niy / 2, niz / 2, fineX, fineY, fineZ, 2.0f);
        System.out.println("{\"source\":\"VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + registerShapesPyramid() + registerShapesExactPyramid()\","
            + "\"transform_bits\":" + vector(atlas.getTransform()) + ",\"lambda_bits\":" + bits((Float)get(atlas, "lambda"))
            + ",\"coarse_initial_energy\":" + coarseInitialEnergy + ",\"coarse_nan_trial_energy\":" + coarseNanTrialEnergy
            + ",\"final_energy\":" + finalEnergy + "}");
    }
}
