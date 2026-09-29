import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import gov.nih.mipav.model.structures.ModelImage;
import java.util.BitSet;

/** Runs selected VandyCRUISE AlgorithmTopologyCorrection on one deterministic tiny 3-D phantom. */
public final class TopologyCorrectionReferenceHarness {
    private static final int NX = 9, NY = 9, NZ = 9;

    public static void main(String[] args) {
        ModelImage source = new ModelImage(new int[] {NX, NY, NZ}, "synthetic_topology_input");
        ModelImage destination = new ModelImage(new int[] {NX, NY, NZ}, "synthetic_topology_output");
        float[] input = new float[NX * NY * NZ];
        BitSet paint = new BitSet(input.length);
        float[] paintValues = new float[input.length];
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    int index = x + NX * y + NX * NY * z;
                    // A deterministic spatial hash provides many distinct levels
                    // for the source fast-marching queue; the corner zero fixes
                    // the normalized range.
                    int hash = (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
                    input[index] = (hash & 0x7fffffff) % 100000 / 100000.0f;
                    if (x == 0 && y == 0 && z == 0) input[index] = 0.0f;
                    if (x == 4 && y == 4 && z == 4) paintValues[index] = 1.0f;
                    if (x == 4 && y == 4 && z == 5) paintValues[index] = 0.5f;
                    if (paintValues[index] > 0.5f) paint.set(index);
                }
            }
        }
        source.setData(input);

        // Same selected 3-D scalar-image call as the Dura adapter's correction:
        // 6/18, paint_mask, 0.5 threshold applied above, [0,1], min/max scaling,
        // 0.05 normalized minimum distance, and downward object->background.
        AlgorithmTopologyCorrection algorithm = new AlgorithmTopologyCorrection(
                destination, source, "scalar_image", 0.05f, 3, 6, 18,
                "paint_mask", paint, 0.0f, 1.0f, true, "object->background");
        algorithm.runAlgorithm();
        if (!algorithm.isCompleted()) throw new IllegalStateException("Algorithm did not complete");

        // Output one xyz-ordered value per line for direct NumPy reshape.
        float[] output = destination.getData();
        System.out.println("OUTPUT_BEGIN");
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    int index = x + NX * y + NX * NY * z;
                    System.out.println(Float.toString(output[index]));
                }
            }
        }
    }
}
