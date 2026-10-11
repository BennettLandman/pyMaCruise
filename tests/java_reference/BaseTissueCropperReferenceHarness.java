import edu.jhmi.rad.medic.utilities.CropParameters;
import edu.jhmi.rad.medic.utilities.CubicVolumeCropper;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.ViewUserInterface;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Direct fixed-input replay of the selected JIST CubicVolumeCropper calls in
 * CropBasedOnGM: derive parameters from atlas GM with threshold 0 and edge 5,
 * then apply those exact parameters to each companion map.
 *
 * Inputs and outputs are big-endian raw arrays in logical NIfTI RAS
 * (rows, cols, slices) order. No orientation or affine transformation is
 * performed by this harness; geometry is checked separately by Python.
 */
public final class BaseTissueCropperReferenceHarness {
    private static final String[] TARGETS = {
        "atlas_gm", "atlas_wm", "atlas_brain_tissue_mask",
        "toads_gm", "toads_wmfill", "toads_mask", "maseg"
    };
    private static final boolean[] INTEGER = {
        true, true, true, false, false, false, true
    };

    private BaseTissueCropperReferenceHarness() {}

    private static ImageData readImage(Path path, int rows, int cols, int slices,
                                       boolean integer) throws IOException {
        ImageData image = integer
            ? new ImageDataInt(rows, cols, slices)
            : new ImageDataFloat(rows, cols, slices);
        try (DataInputStream input = new DataInputStream(
                new BufferedInputStream(Files.newInputStream(path), 1 << 20))) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    for (int k = 0; k < slices; k++) {
                        if (integer) {
                            image.set(i, j, k, input.readInt());
                        } else {
                            image.set(i, j, k, input.readFloat());
                        }
                    }
                }
            }
            if (input.read() != -1) {
                throw new IOException("Input has trailing bytes: " + path);
            }
        } catch (IOException | RuntimeException exception) {
            image.dispose();
            throw exception;
        }
        return image;
    }

    private static void writeImage(Path path, ImageData image,
                                   boolean integer) throws IOException {
        try (DataOutputStream output = new DataOutputStream(
                new BufferedOutputStream(Files.newOutputStream(path), 1 << 20))) {
            for (int i = 0; i < image.getRows(); i++) {
                for (int j = 0; j < image.getCols(); j++) {
                    for (int k = 0; k < image.getSlices(); k++) {
                        if (integer) {
                            output.writeInt(image.getInt(i, j, k));
                        } else {
                            output.writeFloat(image.getFloat(i, j, k));
                        }
                    }
                }
            }
        }
    }

    private static String metadata(CropParameters parameters, int edge,
                                   int rows, int cols, int slices) {
        return "{\n"
            + "  \"source_class\": \"edu.jhmi.rad.medic.utilities.CubicVolumeCropper\",\n"
            + "  \"source_wrapper\": \"CropBasedOnGM: crop(gm, 0, edge), then crop(target, params)\",\n"
            + "  \"threshold\": 0.0,\n"
            + "  \"edge_width\": " + edge + ",\n"
            + "  \"input_shape\": [" + rows + ", " + cols + ", " + slices + "],\n"
            + "  \"bounds_inclusive\": [[" + parameters.xmin + ", " + parameters.xmax
            + "], [" + parameters.ymin + ", " + parameters.ymax + "], ["
            + parameters.zmin + ", " + parameters.zmax + "]],\n"
            + "  \"cropped_shape\": [" + parameters.getCroppedRows() + ", "
            + parameters.getCroppedCols() + ", " + parameters.getCroppedSlices() + "],\n"
            + "  \"background\": " + parameters.background + ",\n"
            + "  \"array_order\": \"logical NIfTI RAS rows, cols, slices\",\n"
            + "  \"header_geometry_checked_in_java\": false\n"
            + "}\n";
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 6) {
            throw new IllegalArgumentException(
                "usage: BaseTissueCropperReferenceHarness <input-dir> <output-dir> "
                    + "<rows> <cols> <slices> <edge-width>");
        }
        Path inputDirectory = Path.of(args[0]);
        Path outputDirectory = Path.of(args[1]);
        int rows = Integer.parseInt(args[2]);
        int cols = Integer.parseInt(args[3]);
        int slices = Integer.parseInt(args[4]);
        int edge = Integer.parseInt(args[5]);
        if (rows <= 0 || cols <= 0 || slices <= 0 || edge < 0) {
            throw new IllegalArgumentException("invalid dimensions or edge width");
        }
        // ImageDataMipav.dispose() logs through JIST preferences. Initialize
        // MIPAV in headless mode and point JIST at a build-local library so
        // its missing-preference fallback never opens a Swing chooser.
        ViewUserInterface.create();
        Path jistLibrary = Path.of(System.getProperty("user.home"), "jist-library");
        Files.createDirectories(jistLibrary);
        Preferences.setProperty("JIST_Library_Directory", jistLibrary.toString());
        Files.createDirectories(outputDirectory);

        CubicVolumeCropper cropper = new CubicVolumeCropper();
        ImageData gm = readImage(inputDirectory.resolve("atlas_gm.i32be"),
                                 rows, cols, slices, true);
        CropParameters parameters;
        ImageData croppedGm = null;
        try {
            croppedGm = cropper.crop(gm, 0.0, edge);
            parameters = cropper.getLastCropParams();
            if (parameters == null) {
                throw new IllegalStateException("Cropper returned no crop parameters");
            }
            writeImage(outputDirectory.resolve("atlas_gm" + (".i32be")),
                       croppedGm, true);
        } finally {
            if (croppedGm != null) croppedGm.dispose();
            gm.dispose();
        }

        for (int index = 1; index < TARGETS.length; index++) {
            String name = TARGETS[index];
            boolean integer = INTEGER[index];
            String suffix = integer ? ".i32be" : ".f32be";
            ImageData target = readImage(inputDirectory.resolve(name + suffix),
                                         rows, cols, slices, integer);
            ImageData cropped = null;
            try {
                cropped = cropper.crop(target, parameters);
                writeImage(outputDirectory.resolve(name + suffix), cropped, integer);
            } finally {
                if (cropped != null) cropped.dispose();
                target.dispose();
            }
        }

        Files.writeString(outputDirectory.resolve("java-crop.json"),
                          metadata(parameters, edge, rows, cols, slices));
    }
}
