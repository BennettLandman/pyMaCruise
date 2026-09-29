import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.vanderbilt.masi.algorithms.CRUISE.ace.MaAnatomicallyConsistentEnhance;

/** Runs the selected MaCRUISE ACE wrapper and VandyCRUISE stages on phantoms. */
public final class ACEReferenceHarness {
    private static void write(DataOutputStream output, ImageData image) throws Exception {
        output.writeInt(image.getRows());
        output.writeInt(image.getCols());
        output.writeInt(image.getSlices());
        for (int x = 0; x < image.getRows(); x++)
            for (int y = 0; y < image.getCols(); y++)
                for (int z = 0; z < image.getSlices(); z++) output.writeFloat(image.getFloat(x, y, z));
    }

    public static void main(String[] args) throws Exception {
        int rows, cols, slices;
        ImageDataFloat gray, white, atlas;
        try (DataInputStream input = new DataInputStream(new FileInputStream(args[0]))) {
            rows = input.readInt(); cols = input.readInt(); slices = input.readInt();
            gray = new ImageDataFloat(rows, cols, slices);
            white = new ImageDataFloat(rows, cols, slices);
            atlas = new ImageDataFloat(rows, cols, slices);
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++) gray.set(x, y, z, input.readFloat());
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++) white.set(x, y, z, input.readFloat());
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++) atlas.set(x, y, z, input.readFloat());
        }
        MaAnatomicallyConsistentEnhance ace = new MaAnatomicallyConsistentEnhance();
        // The optional fourth argument selects the source wrapper's CSF merge
        // placement: 2 after diffusion (base) or 3 after speed Gaussian (plus).
        int placeToAddCSF = args.length > 3 ? Integer.parseInt(args[3]) : 2;
        ace.solve(white, gray, 126.9, new File(args[2]), atlas, placeToAddCSF, 0, 1.8, true);
        try (DataOutputStream output = new DataOutputStream(new FileOutputStream(args[1]))) {
            write(output, ace.getEnhancedGM());
            write(output, ace.getSkeleton());
            write(output, ace.getThinnedSkeleton());
        }
    }
}
