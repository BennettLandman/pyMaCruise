package edu.jhu.ece.iacl.jist.structures.image;

/** Minimal integer image storage required by the original ThinACE source. */
public class ImageDataInt extends ImageData {
    private int[][][] values;

    public ImageDataInt(int x, int y, int z) {
        values = new int[x][y][z];
        rows = x;
        cols = y;
        slices = z;
    }

    public ImageDataInt(ImageData source) {
        this(source.getRows(), source.getCols(), source.getSlices());
        for (int x = 0; x < rows; x++)
            for (int y = 0; y < cols; y++)
                for (int z = 0; z < slices; z++)
                    values[x][y][z] = source.getInt(x, y, z);
        header = source.getHeader();
        name = source.getName();
    }

    public int[][][] toArray3d() { return values; }
    public ImageDataInt mimic() { return new ImageDataInt(rows, cols, slices); }
    public float getFloat(int x, int y, int z) { return values[x][y][z]; }
    public int getInt(int x, int y, int z) { return values[x][y][z]; }
    public boolean getBoolean(int x, int y, int z) { return values[x][y][z] != 0; }
    public byte getByte(int x, int y, int z) { return (byte) values[x][y][z]; }
    public byte getUByte(int x, int y, int z) { return (byte) values[x][y][z]; }
    public double getDouble(int x, int y, int z) { return values[x][y][z]; }
    public void set(int x, int y, int z, int value) { values[x][y][z] = value; }
    public void printNonZeroString(int limit) {}
}
