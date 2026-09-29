import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import gov.nih.mipav.model.structures.ModelImage;
import java.util.BitSet;

/** Runs selected VandyCRUISE topology correction on deterministic small phantoms. */
public final class TopologyCorrectionReferenceHarness {
    private static final class Fixture {
        final int nx, ny, nz;
        final float[] input, paintValues;

        Fixture(int x, int y, int z) {
            nx = x; ny = y; nz = z;
            input = new float[nx * ny * nz];
            paintValues = new float[input.length];
        }

        int index(int x, int y, int z) { return x + nx * y + nx * ny * z; }
        void paint(int x, int y, int z, float value) { paintValues[index(x, y, z)] = value; }
    }

    public static void main(String[] args) {
        String fixtureName = args.length > 0 ? args[0] : "spatial_hash";
        Fixture fixture = makeFixture(fixtureName);
        ModelImage source = new ModelImage(
                new int[] {fixture.nx, fixture.ny, fixture.nz}, "synthetic_topology_input");
        ModelImage destination = new ModelImage(
                new int[] {fixture.nx, fixture.ny, fixture.nz}, "synthetic_topology_output");
        BitSet paint = new BitSet(fixture.input.length);
        for (int index = 0; index < fixture.paintValues.length; index++) {
            if (fixture.paintValues[index] > 0.5f) paint.set(index);
        }
        source.setData(fixture.input);

        // Same selected 3-D scalar-image call as the Dura adapter: 6/18,
        // paint_mask, 0.5 threshold, [0,1], min/max scaling, 0.05 normalized
        // minimum distance, and downward object->background propagation.
        AlgorithmTopologyCorrection algorithm = new AlgorithmTopologyCorrection(
                destination, source, "scalar_image", 0.05f, 3, 6, 18,
                "paint_mask", paint, 0.0f, 1.0f, true, "object->background");
        algorithm.runAlgorithm();
        if (!algorithm.isCompleted()) throw new IllegalStateException("Algorithm did not complete");

        writeSection("INPUT_BEGIN", fixture.input);
        writeSection("PAINT_BEGIN", fixture.paintValues);
        writeSection("OUTPUT_BEGIN", destination.getData());
    }

    private static Fixture makeFixture(String name) {
        switch (name) {
            case "spatial_hash": return spatialHash();
            case "plateau": return plateau();
            case "boundary": return boundary();
            case "narrow_branch": return narrowBranch();
            case "disconnected_seeds": return disconnectedSeeds();
            default: throw new IllegalArgumentException("Unknown fixture: " + name);
        }
    }

    private static Fixture spatialHash() {
        Fixture f = new Fixture(9, 9, 9);
        for (int x = 0; x < f.nx; x++) for (int y = 0; y < f.ny; y++) for (int z = 0; z < f.nz; z++) {
            int hash = (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
            f.input[f.index(x,y,z)] = (hash & 0x7fffffff) % 100000 / 100000.0f;
        }
        f.input[f.index(0,0,0)] = 0.0f;
        f.paint(4, 4, 4, 1.0f);
        f.paint(4, 4, 5, 0.5f);
        return f;
    }

    private static Fixture plateau() {
        Fixture f = new Fixture(9, 9, 9);
        for (int x = 0; x < f.nx; x++) for (int y = 0; y < f.ny; y++) for (int z = 0; z < f.nz; z++) {
            f.input[f.index(x,y,z)] = ((x + 2 * y + 3 * z) % 5) / 4.0f;
        }
        f.input[f.index(0,0,0)] = 0.0f;
        f.input[f.index(8,8,8)] = 1.5f;
        f.input[f.index(0,8,8)] = -0.5f;
        f.paint(4, 4, 4, 1.0f);
        f.paint(4, 4, 5, 0.5f);
        return f;
    }

    private static Fixture boundary() {
        Fixture f = spatialHash();
        for (int i = 0; i < f.paintValues.length; i++) f.paintValues[i] = 0.0f;
        f.paint(0, 4, 4, 1.0f);
        f.paint(0, 4, 5, 0.5f);
        return f;
    }

    private static Fixture narrowBranch() {
        Fixture f = new Fixture(11, 9, 9);
        for (int x = 0; x < f.nx; x++) for (int y = 0; y < f.ny; y++) for (int z = 0; z < f.nz; z++) {
            int dx = x - 5, dy = y - 4, dz = z - 4;
            f.input[f.index(x,y,z)] = Math.min(1.0f, (dx*dx + dy*dy + dz*dz) / 80.0f);
        }
        f.input[f.index(0,0,0)] = 0.0f;
        f.paint(3, 4, 4, 1.0f);
        for (int x = 4; x <= 8; x++) f.paint(x, 4, 4, 1.0f);
        f.paint(8, 4, 5, 0.5f);
        return f;
    }

    private static Fixture disconnectedSeeds() {
        Fixture f = spatialHash();
        for (int i = 0; i < f.paintValues.length; i++) f.paintValues[i] = 0.0f;
        f.paint(2, 2, 2, 1.0f);
        f.paint(6, 6, 6, 1.0f);
        f.paint(2, 2, 3, 0.5f);
        return f;
    }

    private static void writeSection(String heading, float[] values) {
        System.out.println(heading);
        for (float value : values) System.out.println(Float.toString(value));
    }
}
