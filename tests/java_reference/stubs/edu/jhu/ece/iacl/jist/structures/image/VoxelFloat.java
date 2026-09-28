package edu.jhu.ece.iacl.jist.structures.image;
public class VoxelFloat implements Comparable<VoxelFloat> {
    private final float value;
    public VoxelFloat(float value) { this.value = value; }
    public float getFloat() { return value; }
    @Override public int compareTo(VoxelFloat other) { return Float.compare(value, other.value); }
}
