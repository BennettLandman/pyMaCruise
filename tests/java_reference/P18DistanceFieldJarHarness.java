import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;

/**
 * Runs the original DistanceField class from the supplied JIST-Test jar on
 * deterministic small signed fields. This harness intentionally does not
 * compile or call the DistanceFieldyh surrogate.
 */
public final class P18DistanceFieldJarHarness {
    private static final class Phantom {
        final String name;
        final float[][][] input;
        final double maxDistance;

        Phantom(String name, float[][][] input, double maxDistance) {
            this.name = name;
            this.input = input;
            this.maxDistance = maxDistance;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            throw new IllegalArgumentException(
                    "usage: P18DistanceFieldJarHarness ORIGINAL_JIST_JAR OUTPUT_JSON");
        }
        File jar = new File(args[0]).getCanonicalFile();
        File output = new File(args[1]).getCanonicalFile();
        File parent = output.getParentFile();
        if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
            throw new IllegalStateException("cannot create output directory " + parent);
        }

        List<Phantom> cases = new ArrayList<Phantom>();
        cases.add(zeroPlane());
        cases.add(obliquePlane());
        cases.add(boundaryCrossing());
        cases.add(finiteTruncation());

        // Avoid AbstractCalculation.markCompleted(), which tries to initialize
        // MIPAV's UI singleton; solve() itself and all of its numerical calls
        // still dispatch into the original JAR implementation.
        DistanceField solver = new DistanceField() {
            @Override
            public synchronized void markCompleted() {
                // This standalone oracle does not need JIST progress/UI state.
            }
        };
        URL loadedFrom = DistanceField.class.getProtectionDomain().getCodeSource().getLocation();
        String classResourceSha = sha256(DistanceField.class.getResourceAsStream("DistanceField.class"));
        String jarSha = sha256(new FileInputStream(jar));

        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(output), StandardCharsets.UTF_8))) {
            writer.write("{\n  \"reference\": {\n");
            field(writer, "jar_path", jar.getPath(), true, 4);
            field(writer, "jar_sha256", jarSha, true, 4);
            field(writer, "class_name", "edu.jhu.ece.iacl.algorithms.volume.DistanceField", true, 4);
            field(writer, "class_revision", DistanceField.getVersion(), true, 4);
            field(writer, "class_sha256", classResourceSha, true, 4);
            field(writer, "class_loaded_from", loadedFrom.toString(), false, 4);
            writer.write("  },\n  \"array_convention\": \"values are nested in Java [x][y][z] order; distances use unit voxel spacing\",\n");
            writer.write("  \"cases\": [\n");
            for (int index = 0; index < cases.size(); index++) {
                Phantom test = cases.get(index);
                ImageDataFloat source = new ImageDataFloat(test.input);
                source.setName(test.name + "_input");
                ImageDataFloat result = solver.solve(source, test.maxDistance);
                writer.write("    {\n");
                field(writer, "name", test.name, true, 6);
                writer.write("      \"shape_xyz\": [" + test.input.length + ", "
                        + test.input[0].length + ", " + test.input[0][0].length + "],\n");
                writer.write("      \"max_distance_voxels\": " + Double.toString(test.maxDistance) + ",\n");
                writer.write("      \"input_xyz\": ");
                writeArray(writer, test.input, 6);
                writer.write(",\n      \"output_xyz\": ");
                writeArray(writer, result.toArray3d(), 6);
                writer.write("\n    }");
                writer.write(index + 1 == cases.size() ? "\n" : ",\n");
            }
            writer.write("  ]\n}\n");
        }
        System.out.println("Reference class revision: " + DistanceField.getVersion());
        System.out.println("Reference class SHA-256: " + classResourceSha);
        System.out.println("Reference JAR SHA-256: " + jarSha);
        System.out.println("Loaded class from: " + loadedFrom);
        System.out.println("Wrote " + cases.size() + " cases to " + output);
    }

    private static Phantom zeroPlane() {
        float[][][] values = new float[5][5][5];
        for (int x = 0; x < 5; x++) {
            for (int y = 0; y < 5; y++) {
                for (int z = 0; z < 5; z++) values[x][y][z] = x - 2.0f;
            }
        }
        return new Phantom("zero_plane", values, 8.0);
    }

    private static Phantom obliquePlane() {
        float[][][] values = new float[6][5][4];
        for (int x = 0; x < 6; x++) {
            for (int y = 0; y < 5; y++) {
                for (int z = 0; z < 4; z++) {
                    values[x][y][z] = (float) (x + 0.7 * y - 0.45 * z - 2.3);
                }
            }
        }
        return new Phantom("oblique_diagonal_plane", values, 8.0);
    }

    private static Phantom boundaryCrossing() {
        float[][][] values = new float[5][4][3];
        for (int x = 0; x < 5; x++) {
            for (int y = 0; y < 4; y++) {
                for (int z = 0; z < 3; z++) values[x][y][z] = x - 0.25f;
            }
        }
        return new Phantom("interface_near_x_boundary", values, 8.0);
    }

    private static Phantom finiteTruncation() {
        float[][][] values = new float[7][7][7];
        for (int x = 0; x < 7; x++) {
            for (int y = 0; y < 7; y++) {
                for (int z = 0; z < 7; z++) {
                    double dx = x - 3.0;
                    double dy = y - 3.0;
                    double dz = z - 3.0;
                    values[x][y][z] = (float) (Math.sqrt(dx * dx + dy * dy + dz * dz) - 1.2);
                }
            }
        }
        return new Phantom("sphere_with_1p5_voxel_truncation", values, 1.5);
    }

    private static void field(BufferedWriter writer, String key, String value,
                              boolean comma, int indent) throws Exception {
        for (int i = 0; i < indent; i++) writer.write(' ');
        writer.write(quote(key) + ": " + quote(value) + (comma ? ",\n" : "\n"));
    }

    private static void writeArray(BufferedWriter writer, float[][][] values,
                                   int indent) throws Exception {
        writer.write("[\n");
        for (int x = 0; x < values.length; x++) {
            for (int i = 0; i < indent + 2; i++) writer.write(' ');
            writer.write("[\n");
            for (int y = 0; y < values[x].length; y++) {
                for (int i = 0; i < indent + 4; i++) writer.write(' ');
                writer.write('[');
                for (int z = 0; z < values[x][y].length; z++) {
                    if (z > 0) writer.write(", ");
                    writer.write(Float.toString(values[x][y][z]));
                }
                writer.write(']');
                writer.write(y + 1 == values[x].length ? "\n" : ",\n");
            }
            for (int i = 0; i < indent + 2; i++) writer.write(' ');
            writer.write(']');
            writer.write(x + 1 == values.length ? "\n" : ",\n");
        }
        for (int i = 0; i < indent; i++) writer.write(' ');
        writer.write(']');
    }

    private static String quote(String value) {
        return "\"" + value.replace("\\", "\\\\").replace("\"", "\\\"") + "\"";
    }

    private static String sha256(InputStream input) throws Exception {
        MessageDigest digest = MessageDigest.getInstance("SHA-256");
        try (InputStream in = new BufferedInputStream(input)) {
            byte[] buffer = new byte[16384];
            int count;
            while ((count = in.read(buffer)) >= 0) digest.update(buffer, 0, count);
        }
        StringBuilder result = new StringBuilder();
        for (byte value : digest.digest()) result.append(String.format("%02x", value & 0xff));
        return result.toString();
    }
}
