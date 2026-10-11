import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume;
import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume.Interpolation;
import edu.jhu.ece.iacl.algorithms.volume.ReorientVolume.Resolution;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import gov.nih.mipav.model.file.FileIO;
import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.ViewUserInterface;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Runs the supplied source-selected JIST ReorientVolume on saved human inputs.
 *
 * Inputs are a reference NIfTI and label/source/output triples. Each
 * output is the Java ModelImage voxel array in MIPAV's x-fastest order, written
 * as big-endian float32 plus a small JSON geometry record. The caller applies
 * the wrapper's reference-header assignment when comparing geometries.
 */
public final class ReorientVolumeHumanTransferReferenceHarness {
    private static ModelImage read(FileIO io, String path) throws Exception {
        ModelImage image = io.readImage(path);
        if (image == null) throw new IllegalStateException("MIPAV could not read " + path);
        return image;
    }

    private static String quote(String value) {
        return "\"" + value.replace("\\", "\\\\").replace("\"", "\\\"") + "\"";
    }

    private static String ints(int[] values) {
        return Arrays.toString(values).replace(" ", "");
    }

    private static String floats(float[] values) {
        return Arrays.toString(values).replace(" ", "");
    }

    private static String geometry(ModelImage image) {
        FileInfoBase info = image.getFileInfo(0);
        return "{\"extents\":" + ints(image.getExtents())
                + ",\"resolutions\":" + floats(info.getResolutions())
                + ",\"origin\":" + floats(info.getOrigin())
                + ",\"image_orientation\":" + image.getImageOrientation()
                + ",\"axis_orientation\":[" + info.getAxisOrientation(0) + ","
                + info.getAxisOrientation(1) + "," + info.getAxisOrientation(2) + "]}";
    }

    private static void writeResult(ModelImage image, String rawPath) throws Exception {
        int count = image.getDataSize();
        float[] values = new float[count];
        image.exportData(0, count, values);
        try (DataOutputStream out = new DataOutputStream(
                new BufferedOutputStream(new FileOutputStream(rawPath)))) {
            for (float value : values) out.writeFloat(value);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 4 || (args.length - 1) % 3 != 0) {
            throw new IllegalArgumentException(
                    "usage: HARNESS reference.nii.gz label source.nii.gz raw.bin,json ...");
        }
        // FileIO expects the MIPAV singleton, but headless create() does not
        // initialize Swing frames. The launcher redirects user.home to build/.
        ViewUserInterface.create();
        // ImageDataMipav logs through JIST preferences. Give JIST a build-local
        // library directory so the headless logger never opens its chooser.
        java.io.File jistLibrary = new java.io.File(System.getProperty("user.home"), "jist-library");
        jistLibrary.mkdirs();
        Preferences.setProperty("JIST_Library_Directory", jistLibrary.getAbsolutePath());
        FileIO io = new FileIO();
        io.setQuiet(true);
        io.setSuppressProgressBar(true);
        ModelImage referenceModel = read(io, args[0]);
        ImageDataMipav reference = new ImageDataMipav(referenceModel);
        String referenceGeometry = geometry(referenceModel);
        for (int i = 1; i < args.length; i += 3) {
            String label = args[i];
            String sourcePath = args[i + 1];
            String[] outputs = args[i + 2].split(",", -1);
            if (outputs.length != 2) throw new IllegalArgumentException("Expected raw,json output pair");
            ModelImage sourceModel = read(io, sourcePath);
            ImageDataMipav source = new ImageDataMipav(sourceModel);
            String sourceGeometry = geometry(sourceModel);
            ImageDataMipav transferred = ReorientVolume.solve(
                    source, reference, Resolution.UNCHANGED, Interpolation.NEAREST_NEIGHBOR);
            // ToadsCruiseSpaceToOriginalSpace assigns the reference header after solve.
            transferred.setHeader(reference.getHeader());
            ModelImage outputModel = transferred.getModelImageCopy();
            writeResult(outputModel, outputs[0]);
            try (PrintWriter json = new PrintWriter(new FileWriter(outputs[1]))) {
                json.println("{\"label\":" + quote(label)
                        + ",\"source\":" + quote(sourcePath)
                        + ",\"reference\":" + quote(args[0])
                        + ",\"resolution\":\"UNCHANGED\""
                        + ",\"interpolation\":\"NEAREST_NEIGHBOR\""
                        + ",\"source_model_geometry\":" + sourceGeometry
                        + ",\"reference_model_geometry\":" + referenceGeometry
                        + ",\"output_model_geometry_before_reference_header\":" + geometry(transferred.getModelImageDirect())
                        + ",\"output_model_geometry_after_reference_header\":" + geometry(outputModel)
                        + "}");
            }
            outputModel.disposeLocal();
            sourceModel.disposeLocal();
        }
        referenceModel.disposeLocal();
    }
}
