package edu.jhu.ece.iacl.jist.structures.image;

public class ImageDataFloat extends ImageData {
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
    }

    public ImageDataFloat(ImageData image) {
        if (image.vector != null) {
            this.vector = copy(image.vector);
            this.rows = image.rows;
            this.cols = image.cols;
            this.slices = image.slices;
        } else {
            this.scalar = copy(image.scalar);
            this.rows = image.rows;
            this.cols = image.cols;
            this.slices = image.slices;
        }
        header = image.header;
        name = image.name;
    }

    public ImageDataFloat(int x, int y, int z) {
        this(new float[x][y][z]);
    }

    public ImageDataFloat(int x, int y, int z, int components) {
        this(new float[x][y][z][components]);
    }

    public float[][][] toArray3d() {
        return scalar;
    }

    public float[][][][] toArray4d() {
        return vector;
    }

    public void set(int x, int y, int z, float value) {
        scalar[x][y][z] = value;
    }

    private static float[][][] copy(float[][][] values) {
        float[][][] output = new float[values.length][values[0].length][values[0][0].length];
        for (int x = 0; x < values.length; x++)
            for (int y = 0; y < values[0].length; y++)
                System.arraycopy(values[x][y], 0, output[x][y], 0, values[0][0].length);
        return output;
    }

    private static float[][][][] copy(float[][][][] values) {
        int components = values[0][0][0].length;
        float[][][][] output = new float[values.length][values[0].length][values[0][0].length][components];
        for (int x = 0; x < values.length; x++)
            for (int y = 0; y < values[0].length; y++)
                for (int z = 0; z < values[0][0].length; z++)
                    System.arraycopy(values[x][y][z], 0, output[x][y][z], 0, components);
        return output;
    }
}
