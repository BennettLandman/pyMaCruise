import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Human-input replay of the selected
 * LesionToadSegmentation.generateDuraRemovalOutputs method. The copied
 * accumulation below follows VandyCRUISE source SHA-256
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * This runs the selected method body on held-fixed arrays; it is not the
 * complete MIPAV/JIST LesionTOADS wrapper.
 *
 * Input binary: big-endian nx, ny, nz, class count; class-major float32
 * memberships; x/y/z-order uint8 segmentation; float32 lesion membership;
 * float32 black-hole membership. Output binary: four channel-major float32
 * products in source order (CSF, GM, WM-fill, WM-mask).
 */
public final class LesionToadsHumanDuraOutputsReferenceHarness {
    private static final String[] NAMES = {
        "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM",
        "Cerebrum-GM", "Caudate", "Thalamus", "Putamen", "Brainstem",
        "Cerebellum-WM", "OuterVentricular-WM", "InterVentricular-WM",
        "Frontal-WM"
    };

    private static boolean csf(String name) {
        return name.equals("SulcalCSF") || name.equals("Sulcal-CSF");
    }

    private static boolean gm(String name) {
        return name.equals("CorticalGM") || name.equals("Cerebrum-GM")
            || name.equals("CerebralGM") || name.equals("Cerebellum-GM")
            || name.equals("CerebellarGM");
    }

    private static boolean excludedFromWm(String name) {
        return csf(name) || name.equals("CorticalGM")
            || name.equals("Cerebrum-GM") || name.equals("CerebralGM")
            || name.equals("Cerebellum-GM") || name.equals("CerebellarGM")
            || name.equals("Background");
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("usage: INPUT OUTPUT");
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(args[0]), 1 << 20))) {
            int nx = in.readInt(), ny = in.readInt(), nz = in.readInt(), classes = in.readInt();
            if (classes != NAMES.length) throw new IllegalArgumentException("expected 13 classes");
            int n = Math.multiplyExact(Math.multiplyExact(nx, ny), nz);
            float[][] memberships = new float[classes][n];
            for (int k = 0; k < classes; k++)
                for (int v = 0; v < n; v++) memberships[k][v] = in.readFloat();
            byte[] labels = new byte[n]; in.readFully(labels);
            float[] lesion = new float[n], blackHole = new float[n];
            for (int v = 0; v < n; v++) lesion[v] = in.readFloat();
            for (int v = 0; v < n; v++) blackHole[v] = in.readFloat();
            if (in.read() != -1) throw new IllegalArgumentException("unexpected trailing input");

            float[][] output = new float[4][n];
            for (int k = 0; k < classes; k++) if (csf(NAMES[k]))
                System.arraycopy(memberships[k], 0, output[0], 0, n);
            for (int v = 0; v < n; v++) output[2][v] = lesion[v] + blackHole[v];
            for (int k = 0; k < classes; k++) if (!excludedFromWm(NAMES[k])) {
                for (int v = 0; v < n; v++) {
                    output[2][v] += memberships[k][v];
                    if ((labels[v] & 0xff) == k) output[3][v] = 1.0f;
                }
            }
            for (int k = 0; k < classes; k++) if (gm(NAMES[k])) {
                for (int v = 0; v < n; v++) {
                    output[1][v] += memberships[k][v];
                    String selectedName = NAMES[labels[v] & 0xff];
                    if (selectedName.equals("Cerebellum-GM")
                            || selectedName.equals("CerebellarGM)")) output[3][v] = 1.0f;
                }
            }
            try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(
                        new FileOutputStream(args[1]), 1 << 20))) {
                for (float[] channel : output) for (float value : channel) out.writeFloat(value);
            }
        }
    }
}
