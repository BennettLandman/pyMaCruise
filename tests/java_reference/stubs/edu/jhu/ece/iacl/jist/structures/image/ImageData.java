package edu.jhu.ece.iacl.jist.structures.image;

/** Minimal scalar/vector image API shared by the TGDM, ThinACE and GVF harnesses. */
public class ImageData {
    protected float[][][] scalar;
    protected float[][][][] vector;
    protected int rows;
    protected int cols;
    protected int slices;
    protected ImageHeader header = new ImageHeader();
    protected String name = "synthetic";

    public int getRows() { return rows; }
    public int getCols() { return cols; }
    public int getSlices() { return slices; }
    public int getComponents() { return vector == null ? 1 : vector[0][0][0].length; }
    public float getFloat(int x, int y, int z) { return scalar[x][y][z]; }
    public float getFloat(int x, int y, int z, int component) {
        return vector == null ? scalar[x][y][z] : vector[x][y][z][component];
    }
    public int getInt(int x, int y, int z) { return (int) getFloat(x, y, z); }
    public byte getByte(int x, int y, int z) { return (byte) getInt(x, y, z); }
    public byte getUByte(int x, int y, int z) { return (byte) getInt(x, y, z); }
    public boolean getBoolean(int x, int y, int z) { return getInt(x, y, z) != 0; }
    public double getDouble(int x, int y, int z) { return getFloat(x, y, z); }
    public void set(int x, int y, int z, float value) { scalar[x][y][z] = value; }
    public void set(int x, int y, int z, double value) { scalar[x][y][z] = (float) value; }
    public void set(int x, int y, int z, int value) { scalar[x][y][z] = value; }
    public void set(int x, int y, int z, int component, float value) { vector[x][y][z][component] = value; }
    public ImageHeader getHeader() { return header; }
    public void setHeader(ImageHeader value) { header = value; }
    public String getName() { return name; }
    public void setName(String value) { name = value; }
}
