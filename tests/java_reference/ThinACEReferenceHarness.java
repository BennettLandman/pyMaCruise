import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import edu.jhu.ece.iacl.algorithms.ace.ThinACE;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;

/** Runs the supplied VandyCRUISE ThinACE source on Python-generated fixtures. */
public final class ThinACEReferenceHarness {
    public static void main(String[] args) throws Exception {
        String inputPath = args[0];
        String outputPath = args[1];
        try (DataInputStream input = new DataInputStream(new FileInputStream(inputPath))) {
            int rows = input.readInt();
            int cols = input.readInt();
            int slices = input.readInt();
            ImageDataInt volume = new ImageDataInt(rows, cols, slices);
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++)
                        volume.set(x, y, z, input.readInt());
            ImageData output = new ThinACE().solve(volume);
            try (DataOutputStream result = new DataOutputStream(new FileOutputStream(outputPath))) {
                result.writeInt(rows);
                result.writeInt(cols);
                result.writeInt(slices);
                for (int x = 0; x < rows; x++)
                    for (int y = 0; y < cols; y++)
                        for (int z = 0; z < slices; z++)
                            result.writeInt(output.getInt(x, y, z));
            }
            if (args.length > 2) {
                ImageDataFloat initialLevelSet = new ImageDataFloat(rows, cols, slices);
                for (int x = 0; x < rows; x++)
                    for (int y = 0; y < cols; y++)
                        for (int z = 0; z < slices; z++)
                            initialLevelSet.set(x, y, z, volume.getInt(x, y, z) > 0 ? 1.0f : -1.0f);
                ImageDataFloat distance = new DistanceField().solve(initialLevelSet, 4.0);
                try (DataOutputStream result = new DataOutputStream(new FileOutputStream(args[2]))) {
                    for (int x = 0; x < rows; x++)
                        for (int y = 0; y < cols; y++)
                            for (int z = 0; z < slices; z++)
                                result.writeFloat(distance.toArray3d()[x][y][z]);
                }
            }
        }
    }
}
