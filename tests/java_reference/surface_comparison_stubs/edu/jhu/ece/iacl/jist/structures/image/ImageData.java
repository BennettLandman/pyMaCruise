package edu.jhu.ece.iacl.jist.structures.image;

/** Minimal scalar-image API needed to execute IsoSurfaceOnGrid.solveOriginal. */
public class ImageData {
    protected int rows;
    protected int cols;
    protected int slices;
    protected String name = "surface-comparison";

    public int getRows() { return rows; }
    public int getCols() { return cols; }
    public int getSlices() { return slices; }
    public String getName() { return name; }
}
