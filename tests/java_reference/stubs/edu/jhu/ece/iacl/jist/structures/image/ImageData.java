package edu.jhu.ece.iacl.jist.structures.image;

/** Minimal scalar/vector image API used by the original TGDM implementation. */
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
    public float getFloat(int x, int y, int z) { return scalar[x][y][z]; }
    public ImageHeader getHeader() { return header; }
    public void setHeader(ImageHeader value) { header = value; }
    public String getName() { return name; }
    public void setName(String value) { name = value; }
}
