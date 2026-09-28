import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import edu.jhu.ece.iacl.algorithms.gvf.GradVecFlowGeneralized;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/** Runs the supplied generalized VandyCRUISE GVF source on a synthetic array. */
public final class GVFReferenceHarness {
    public static void main(String[] args) throws Exception {
        int rows, cols, slices;
        ImageDataFloat input;
        try (DataInputStream stream = new DataInputStream(new FileInputStream(args[0]))) {
            rows = stream.readInt(); cols = stream.readInt(); slices = stream.readInt();
            input = new ImageDataFloat(rows, cols, slices);
            float[][][] values = input.toArray3d();
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++) values[x][y][z] = stream.readFloat();
        }
        float lambda = Float.parseFloat(args[2]);
        int iterations = Integer.parseInt(args[3]);
        ImageDataFloat result = (ImageDataFloat) new GradVecFlowGeneralized().solve(input, lambda, false, iterations);
        try (DataOutputStream stream = new DataOutputStream(new FileOutputStream(args[1]))) {
            stream.writeInt(rows); stream.writeInt(cols); stream.writeInt(slices);
            float[][][][] vectors = result.toArray4d();
            for (int x = 0; x < rows; x++)
                for (int y = 0; y < cols; y++)
                    for (int z = 0; z < slices; z++)
                        for (int c = 0; c < 3; c++) stream.writeFloat(vectors[x][y][z][c]);
        }
    }
}
