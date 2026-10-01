package edu.jhu.ece.iacl.jist.structures.image;

/** Float32 scalar image adapter with the methods used by IsoSurfaceOnGrid. */
public final class ImageDataFloat extends ImageData {
    private final float[][][] data;

    public ImageDataFloat(float[][][] values) {
        rows = values.length;
        cols = values[0].length;
        slices = values[0][0].length;
        data = new float[rows][cols][slices];
        for (int x = 0; x < rows; x++)
            for (int y = 0; y < cols; y++)
                System.arraycopy(values[x][y], 0, data[x][y], 0, slices);
    }

    public ImageDataFloat(ImageData source) {
        ImageDataFloat copy = (ImageDataFloat) source;
        rows = copy.rows;
        cols = copy.cols;
        slices = copy.slices;
        name = copy.name;
        data = new float[rows][cols][slices];
        for (int x = 0; x < rows; x++)
            for (int y = 0; y < cols; y++)
                System.arraycopy(copy.data[x][y], 0, data[x][y], 0, slices);
    }

    public float[][][] toArray3d() { return data; }
}
