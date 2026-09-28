package edu.jhu.ece.iacl.jist.io;

import java.io.File;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;

/** Intermediate-image writing is disabled in the synthetic ACE reference run. */
public final class ImageDataReaderWriter {
    private static final ImageDataReaderWriter INSTANCE = new ImageDataReaderWriter();
    private ImageDataReaderWriter() {}
    public static ImageDataReaderWriter getInstance() { return INSTANCE; }
    public void write(ImageData image, File directory) {
        try {
            directory.mkdirs();
            File output = new File(directory, image.getName() + ".bin");
            try (DataOutputStream stream = new DataOutputStream(new FileOutputStream(output))) {
                stream.writeInt(image.getRows());
                stream.writeInt(image.getCols());
                stream.writeInt(image.getSlices());
                for (int x = 0; x < image.getRows(); x++)
                    for (int y = 0; y < image.getCols(); y++)
                        for (int z = 0; z < image.getSlices(); z++) stream.writeFloat(image.getFloat(x, y, z));
            }
        } catch (Exception error) {
            throw new RuntimeException("Failed to capture synthetic ACE intermediate", error);
        }
    }
}
