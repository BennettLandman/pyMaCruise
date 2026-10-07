import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

/** Runs the selected Dura wrapper sequence on held-fixed workflow inputs. */
public final class DuraBaseInputsReferenceHarness {
    private static final int NX = Integer.getInteger("dura.nx", 172);
    private static final int NY = Integer.getInteger("dura.ny", 216);
    private static final int NZ = Integer.getInteger("dura.nz", 164);

    private static float[][][] readVolume(DataInputStream in) throws Exception {
        float[][][] values = new float[NX][NY][NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            values[x][y][z] = in.readFloat();
        return values;
    }

    private static void writeVolume(Path path, float[][][] values) throws Exception {
        try (DataOutputStream out = new DataOutputStream(new GZIPOutputStream(
                new FileOutputStream(path.toFile()), 1 << 20))) {
            out.writeInt(NX); out.writeInt(NY); out.writeInt(NZ);
            for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
                out.writeFloat(values[x][y][z]);
        }
    }

    private static void writeMask(Path path, byte[][][] values) throws Exception {
        try (DataOutputStream out = new DataOutputStream(new GZIPOutputStream(
                new FileOutputStream(path.toFile()), 1 << 20))) {
            out.writeInt(NX); out.writeInt(NY); out.writeInt(NZ);
            for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
                out.writeByte(values[x][y][z]);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("usage: INPUT_FLOAT32_BE OUTPUT_DIR");
        Path input = Path.of(args[0]);
        Path output = Path.of(args[1]);
        Files.createDirectories(output);
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(input.toFile()), 1 << 20))) {
            int nx = in.readInt(), ny = in.readInt(), nz = in.readInt();
            if (nx != NX || ny != NY || nz != NZ) throw new IllegalArgumentException("input dimensions disagree with JVM properties");
            float[][][] t1 = readVolume(in), wm = readVolume(in), gm = readVolume(in), mask = readVolume(in);
            if (in.read() != -1) throw new IllegalArgumentException("unexpected trailing input bytes");
            DuraAttachedNoDuraReferenceHarness.CaseResult result =
                DuraAttachedNoDuraReferenceHarness.runFromInputs("base_toads_first", t1, wm, gm, mask);
            writeVolume(output.resolve("wm_tc.raw.gz"), result.wmTc);
            writeVolume(output.resolve("union.raw.gz"), result.union);
            writeVolume(output.resolve("gm_tc.raw.gz"), result.gmTc);
            writeVolume(output.resolve("inner_distance.raw.gz"), result.innerDistance);
            writeVolume(output.resolve("outer_distance.raw.gz"), result.outerDistance);
            writeVolume(output.resolve("level_set.raw.gz"), result.levelSet);
            writeVolume(output.resolve("stripped_t1.raw.gz"), result.stripped);
            writeVolume(output.resolve("corrected_wm.raw.gz"), result.wmOut);
            writeVolume(output.resolve("corrected_gm.raw.gz"), result.gmOut);
            writeVolume(output.resolve("brain_mask_level_set.raw.gz"), result.maskLevelSet);
            writeMask(output.resolve("brain_mask.raw.gz"), result.mask);
            for (int i=0;i<result.evolution.passLevelSets.size();i++)
                writeVolume(output.resolve(String.format("tgdm_pass_%02d.raw.gz", i+1)), result.evolution.passLevelSets.get(i));
        }
    }
}
