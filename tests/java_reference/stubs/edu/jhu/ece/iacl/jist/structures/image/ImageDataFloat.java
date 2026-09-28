package edu.jhu.ece.iacl.jist.structures.image;

/** Small in-memory subset of JIST ImageDataFloat shared by numerical harnesses. */
public class ImageDataFloat extends ImageData implements Cloneable {
    private int components = 1;

    public ImageDataFloat(float[][][] data) {
        scalar = data;
        rows = data.length;
        cols = data[0].length;
        slices = data[0][0].length;
    }

    public ImageDataFloat(float[][][][] data) {
        vector = data;
        rows = data.length;
        cols = data[0].length;
        slices = data[0][0].length;
        components = data[0][0][0].length;
    }

    public ImageDataFloat(int rows, int cols, int slices) {
        this(new float[rows][cols][slices]);
    }

    public ImageDataFloat(int rows, int cols, int slices, int components) {
        this(new float[rows][cols][slices][components]);
    }

    public ImageDataFloat(ImageData source) {
        rows = source.getRows();
        cols = source.getCols();
        slices = source.getSlices();
        components = source.getComponents();
        if (source.vector == null) scalar = copy(source.scalar);
        else vector = copy(source.vector);
        header = source.getHeader();
        name = source.getName();
    }

    public float[][][] toArray3d() {
        if (scalar == null) throw new IllegalStateException("not a scalar image");
        return scalar;
    }

    public float[][][][] toArray4d() {
        if (vector == null) throw new IllegalStateException("not a vector image");
        return vector;
    }

    public float getFloat(int x, int y, int z, int component) {
        return vector == null ? scalar[x][y][z] : vector[x][y][z][component];
    }

    public void set(int x, int y, int z, float value) { scalar[x][y][z] = value; }
    public void set(int x, int y, int z, int component, float value) { vector[x][y][z][component] = value; }
    @Override public int getComponents() { return components; }

    public ImageDataFloat mimic() {
        return components == 1
            ? new ImageDataFloat(rows, cols, slices)
            : new ImageDataFloat(rows, cols, slices, components);
    }

    public void normalize() {
        if (vector == null) throw new IllegalStateException("only vector images can be normalized");
        for (int x = 0; x < rows; x++)
            for (int y = 0; y < cols; y++)
                for (int z = 0; z < slices; z++) {
                    float length = 0;
                    for (int c = 0; c < components; c++) length += vector[x][y][z][c] * vector[x][y][z][c];
                    if (length > 0) {
                        length = (float) Math.sqrt(length);
                        for (int c = 0; c < components; c++) vector[x][y][z][c] /= length;
                    }
                }
    }

    @Override public ImageDataFloat clone() { return new ImageDataFloat(this); }

    private static float[][][] copy(float[][][] source) {
        float[][][] out = new float[source.length][source[0].length][source[0][0].length];
        for (int x = 0; x < source.length; x++)
            for (int y = 0; y < source[0].length; y++)
                System.arraycopy(source[x][y], 0, out[x][y], 0, source[0][0].length);
        return out;
    }

    private static float[][][][] copy(float[][][][] source) {
        int n = source[0][0][0].length;
        float[][][][] out = new float[source.length][source[0].length][source[0][0].length][n];
        for (int x = 0; x < source.length; x++)
            for (int y = 0; y < source[0].length; y++)
                for (int z = 0; z < source[0][0].length; z++)
                    System.arraycopy(source[x][y][z], 0, out[x][y][z], 0, n);
        return out;
    }
}
