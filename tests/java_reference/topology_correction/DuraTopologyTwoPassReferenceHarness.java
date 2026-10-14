import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import gov.nih.mipav.model.structures.ModelImage;
import java.util.BitSet;

/** Runs the two selected topology calls and float32 union from RemoveDura09. */
public final class DuraTopologyTwoPassReferenceHarness {
    private static final int NX = 9, NY = 9, NZ = 9;

    public static void main(String[] args) {
        int size = NX * NY * NZ;
        float[] wm = new float[size];
        float[] gm = new float[size];
        float[] paintValues = new float[size];
        BitSet paint = new BitSet(size);
        for (int x = 0; x < NX; x++) for (int y = 0; y < NY; y++) for (int z = 0; z < NZ; z++) {
            int index = x + NX * y + NX * NY * z;
            int hash = (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
            wm[index] = (hash & 0x7fffffff) % 100000 / 100000.0f;
            if (x == 0 && y == 0 && z == 0) wm[index] = 0.0f;
            gm[index] = 0.2f + ((x * 3 + y * 5 + z * 7) % 23) / 100.0f;
            if (x == 4 && y == 4 && z == 4) paintValues[index] = 1.0f;
            if (x == 4 && y == 4 && z == 5) paintValues[index] = 0.5f;
            if (paintValues[index] > 0.5f) paint.set(index);
        }

        float[] wmTC = correct(wm, paint, "wm", 0.05f);
        float[] union = new float[size];
        for (int i = 0; i < size; i++) union[i] = wmTC[i] + gm[i];
        float[] gmTC = correct(union, paint, "gm", 0.02f);

        writeSection("WM_INPUT_BEGIN", wm);
        writeSection("PAINT_BEGIN", paintValues);
        writeSection("RAW_GM_BEGIN", gm);
        writeSection("WM_TC_BEGIN", wmTC);
        writeSection("UNION_BEGIN", union);
        writeSection("GM_TC_BEGIN", gmTC);
    }

    private static float[] correct(float[] values, BitSet paint, String name, float regularize) {
        ModelImage source = new ModelImage(new int[] {NX, NY, NZ}, name + "_input");
        ModelImage destination = new ModelImage(new int[] {NX, NY, NZ}, name + "_output");
        source.setData(values);
        AlgorithmTopologyCorrection algorithm = new AlgorithmTopologyCorrection(
                destination, source, "scalar_image", regularize, 3, 6, 18,
                "paint_mask", paint, 0.0f, 1.0f, true, "object->background");
        algorithm.runAlgorithm();
        if (!algorithm.isCompleted()) throw new IllegalStateException("Algorithm did not complete: " + name);
        return destination.getData();
    }

    private static void writeSection(String heading, float[] values) {
        System.out.println(heading);
        for (float value : values) System.out.println(Float.toString(value));
    }
}
