import edu.jhmi.rad.medic.algorithms.AlgorithmMultipleToads;
import edu.jhmi.rad.medic.methods.MultipleToadDeformableAtlas;
import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.ModelStorageBase;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/**
 * Runs the selected V3 TOADS wrapper and exports its classification and four
 * DuraRemoval_input products in ModelImage x-fastest order. The selected
 * AlgorithmMultipleToads and MultipleToadSegmentation source are compiled
 * unchanged; output conversion is left to a separate Python comparator.
 */
public final class ToadsWorkflowDuraOutputsHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 14) {
            throw new IllegalArgumentException(
                    "input.raw atlas.txt nx ny nz sx sy sz maxIterations classification.raw csf.raw gm.raw wmfill.raw wmmask.raw");
        }
        int nx = Integer.parseInt(args[2]);
        int ny = Integer.parseInt(args[3]);
        int nz = Integer.parseInt(args[4]);
        int maxIterations = Integer.parseInt(args[8]);
        byte[] raw = Files.readAllBytes(Path.of(args[0]));
        if (raw.length != 4 * nx * ny * nz) {
            throw new IllegalArgumentException("wrong raw T1 length");
        }
        ByteBuffer input = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[nx * ny * nz];
        for (int i = 0; i < values.length; i++) values[i] = input.getFloat();
        ModelImage image = new ModelImage(
                ModelStorageBase.FLOAT, new int[]{nx, ny, nz}, "synthetic T1");
        image.importData(0, values, true);
        FileInfoBase info = image.getFileInfo(0);
        info.setResolutions(new float[]{
                Float.parseFloat(args[5]), Float.parseFloat(args[6]), Float.parseFloat(args[7])});
        info.setUnitsOfMeasure(new int[]{
                FileInfoBase.MILLIMETERS, FileInfoBase.MILLIMETERS, FileInfoBase.MILLIMETERS});
        info.setAxisOrientation(new int[]{
                FileInfoBase.ORI_L2R_TYPE, FileInfoBase.ORI_P2A_TYPE, FileInfoBase.ORI_I2S_TYPE});
        info.setImageOrientation(FileInfoBase.AXIAL);
        image.setImageOrientation(FileInfoBase.AXIAL);

        MultipleToadDeformableAtlas atlas = new MultipleToadDeformableAtlas(args[1]);
        AlgorithmMultipleToads algorithm = new AlgorithmMultipleToads(
                new ModelImage[]{image}, 1, new String[]{"T1_MPRAGE"}, args[1], atlas,
                "DuraRemoval_input", 0.2f, maxIterations, 1.0e-4f, 3.0f, 0.5f,
                2.0f, 0.2f, "rigid", 4, 50, 5, 0, true, "Chebyshev", 3, 30.0f, "18/6");
        algorithm.runAlgorithm();
        ModelImage[] outputs = algorithm.getResultImages();
        if (!algorithm.isCompleted() || outputs == null || outputs.length != 5) {
            throw new IllegalStateException("Java TOADS wrapper did not produce five Dura inputs");
        }
        for (int output = 0; output < outputs.length; output++) {
            if (outputs[output] == null || outputs[output].getDataSize() != nx * ny * nz) {
                throw new IllegalStateException("invalid TOADS output index " + output);
            }
            ByteBuffer bytes = ByteBuffer.allocate(nx * ny * nz * (output == 0 ? 1 : 4))
                    .order(ByteOrder.LITTLE_ENDIAN);
            if (output == 0) {
                byte[] labels = new byte[nx * ny * nz];
                outputs[output].exportData(0, labels.length, labels);
                bytes.put(labels);
            } else {
                float[] field = new float[nx * ny * nz];
                outputs[output].exportData(0, field.length, field);
                for (float value : field) bytes.putFloat(value);
            }
            Files.write(Path.of(args[9 + output]), bytes.array());
        }
        System.out.println("completed=" + algorithm.isCompleted() + " dimensions="
                + Arrays.toString(outputs[0].getExtents()) + " exported_products=" + outputs.length);
        image.disposeLocal();
        for (ModelImage result : outputs) result.disposeLocal();
    }
}
