package gov.nih.mipav.model.structures;

import java.io.IOException;

/** Float-array image shim for the in-memory topology source harness. */
public class ModelImage {
    private final int[] extents;
    private final String name;
    private float[] values;
    public ModelImage(int[] dimensions, String imageName) {
        extents = dimensions.clone(); name = imageName;
        values = new float[dimensions[0] * dimensions[1] * dimensions[2]];
    }
    public int getNDims() { return extents.length; }
    public int[] getExtents() { return extents.clone(); }
    public String getImageName() { return name; }
    public void exportData(int start, int length, float[] target) throws IOException {
        System.arraycopy(values, start, target, 0, length);
    }
    public void importData(int start, float[] source, boolean update) throws IOException {
        System.arraycopy(source, 0, values, start, source.length);
    }
    public float[] getData() { return values.clone(); }
    public void setData(float[] source) { values = source.clone(); }
    public Object getHistoryArea() { return null; }
}
