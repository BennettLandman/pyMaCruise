import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;

/** Exercises the Dura wrapper arithmetic and original JIST distance solver. */
public final class DuraInitialLevelSetsJarHarness {
    private static final int NX = 9, NY = 8, NZ = 7;
    private static final float MEMBERSHIP = 0.4f;

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("usage: OUTPUT_JSON");
        float[][][] wm = new float[NX][NY][NZ];
        float[][][] gm = new float[NX][NY][NZ];
        float[][][] inner = new float[NX][NY][NZ];
        float[][][] outer = new float[NX][NY][NZ];
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++)
            for (int z = 0; z < NZ; z++) {
                wm[x][y][z] = (float)(132.0 + 74.0 * Math.sin(0.41*x + 0.17*y - 0.2*z));
                gm[x][y][z] = (float)(104.0 + 86.0 * Math.cos(0.19*x - 0.33*y + 0.27*z));
                inner[x][y][z] = 0.5f - wm[x][y][z] / 255.0f;
                outer[x][y][z] = MEMBERSHIP - gm[x][y][z] / 255.0f;
            }
        DistanceField solver = new DistanceField() {
            @Override public synchronized void markCompleted() { }
        };
        float[][][] innerDistance = solver.solve(new ImageDataFloat(inner), 15.0).toArray3d();
        float[][][] outerDistance = solver.solve(new ImageDataFloat(outer), 10.0).toArray3d();

        try (BufferedWriter out = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(new File(args[0])), StandardCharsets.UTF_8))) {
            out.write("{\"source\":\"MedicAlgorithmRemoveDura09 arithmetic + original JIST DistanceField\",");
            out.write("\"jar_distancefield_sha256\":\"92294be2db74351446c3e5898fc35b58497f89ccb38772d4d565885f3f0c35b7\",");
            out.write("\"shape_xyz\":["+NX+","+NY+","+NZ+"],\"membership_level\":"+MEMBERSHIP+",");
            field(out, "white_matter", wm); out.write(',');
            field(out, "gray_matter", gm); out.write(',');
            field(out, "inner_level_set", inner); out.write(',');
            field(out, "outer_level_set", outer); out.write(',');
            field(out, "inner_signed_distance", innerDistance); out.write(',');
            field(out, "outer_signed_distance", outerDistance);
            out.write('}');
        }
    }

    private static void field(BufferedWriter out, String name, float[][][] values) throws Exception {
        out.write("\"" + name + "\":[");
        boolean first = true;
        for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) {
            if (!first) out.write(',');
            first = false;
            out.write(Float.toString(values[x][y][z]));
        }
        out.write(']');
    }
}
