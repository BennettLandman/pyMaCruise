import edu.jhmi.rad.medic.methods.DemonToadDeformableAtlas;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.lang.reflect.Field;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

/** Executes the supplied AtlasEM rigid optimizer on caller-provided float32 arrays. */
public final class AtlasEMDemonRigidRegistrationReferenceHarness {
    private AtlasEMDemonRigidRegistrationReferenceHarness() {}

    private static float[] readFloats(Path path, int count) throws IOException {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length != count * Float.BYTES) {
            throw new IllegalArgumentException(path + " has " + bytes.length + " bytes; expected " + (count * Float.BYTES));
        }
        FloatBufferView view = new FloatBufferView(bytes);
        float[] values = new float[count];
        for (int i = 0; i < count; i++) values[i] = view.get(i);
        return values;
    }

    private static final class FloatBufferView {
        private final ByteBuffer buffer;
        FloatBufferView(byte[] bytes) { buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN); }
        float get(int index) { return buffer.getFloat(index * Float.BYTES); }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 7) throw new IllegalArgumentException("nx ny nz iterations levels single|pyramid arrays-directory");
        int nx = Integer.parseInt(args[0]);
        int ny = Integer.parseInt(args[1]);
        int nz = Integer.parseInt(args[2]);
        int iterations = Integer.parseInt(args[3]);
        int levels = Integer.parseInt(args[4]);
        String mode = args[5];
        Path directory = Path.of(args[6]);
        int voxels = nx * ny * nz;
        byte[] shapeBytes = Files.readAllBytes(directory.resolve("shapes.f32le"));
        if (shapeBytes.length % (voxels * Float.BYTES) != 0) {
            throw new IllegalArgumentException("shape payload is not an integral class count");
        }
        int classes = shapeBytes.length / (voxels * Float.BYTES);
        float[] shapeValues = readFloats(directory.resolve("shapes.f32le"), classes * voxels);
        float[] membershipValues = readFloats(directory.resolve("memberships.f32le"), classes * voxels);

        float[][][][] shapes = new float[classes][nx][ny][nz];
        float[][][][] memberships = new float[nx][ny][nz][classes];
        int offset = 0;
        for (int k = 0; k < classes; k++)
            for (int x = 0; x < nx; x++)
                for (int y = 0; y < ny; y++)
                    for (int z = 0; z < nz; z++) shapes[k][x][y][z] = shapeValues[offset++];
        offset = 0;
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++)
                    for (int k = 0; k < classes; k++) memberships[x][y][z][k] = membershipValues[offset++];

        PrintStream standardOutput = System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));
        DemonToadDeformableAtlas atlas = new DemonToadDeformableAtlas(classes);
        atlas.setShapes(shapes, new int[] {nx, ny, nz}, new float[] {0.8f, 1.0f, 1.25f});
        atlas.setImageInfo(nx, ny, nz, 1.0f, 1.2f, 1.5f, 0, 0, 1, 2);
        for (int k = 0; k < classes; k++) {
            atlas.setName(k, k == 0 ? "background" : "class-" + k);
            atlas.setLabel(k, (byte) k);
        }
        boolean[] registered = new boolean[classes];
        for (int k = 1; k < classes; k++) registered[k] = true;
        Path registeredPath = directory.resolve("registered.u8");
        if (Files.exists(registeredPath)) {
            byte[] flags = Files.readAllBytes(registeredPath);
            if (flags.length != classes) throw new IllegalArgumentException("registered.u8 has wrong class count");
            for (int k = 0; k < classes; k++) registered[k] = flags[k] != 0;
        }
        Field registeredField = DemonToadDeformableAtlas.class.getDeclaredField("registeredShape");
        registeredField.setAccessible(true);
        registeredField.set(atlas, registered);

        atlas.initShapeRegistration(memberships, "rigid", iterations, levels);
        if (mode.equals("single")) atlas.registerShapes();
        else if (mode.equals("pyramid")) atlas.registerShapesPyramid();
        else throw new IllegalArgumentException("mode must be single or pyramid");
        float[] transform = atlas.getTransform();
        System.setOut(standardOutput);
        StringBuilder output = new StringBuilder();
        for (int i = 0; i < transform.length; i++) {
            if (i > 0) output.append(',');
            output.append(Float.floatToIntBits(transform[i]));
        }
        System.out.println(output);
    }
}
