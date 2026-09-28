package edu.jhu.ece.iacl.jist.structures.image;

import edu.jhu.ece.iacl.jist.structures.data.Indexable;

public final class VoxelIndexed<T extends VoxelFloat> implements Indexable<VoxelFloat> {
    private T value;
    private int row, column, slice;
    private int chainIndex, index;
    public VoxelIndexed(T value) { this.value = value; }
    public void setRefPosition(int row, int column, int slice) {
        this.row = row; this.column = column; this.slice = slice;
    }
    public int getRow() { return row; }
    public int getColumn() { return column; }
    public int getSlice() { return slice; }
    public float getFloat() { return value.getFloat(); }
    public double getDouble() { return value.getDouble(); }
    public int getInt() { return value.getInt(); }
    public int getChainIndex() { return chainIndex; }
    public void setChainIndex(int value) { chainIndex = value; }
    public void setIndex(int value) { index = value; }
    public int getIndex() { return index; }
    public VoxelFloat getValue() { return value; }
    @SuppressWarnings("unchecked")
    public void setValue(Comparable value) { this.value = (T) value; }
    @SuppressWarnings("unchecked")
    public void set(double value) { this.value = (T) new VoxelFloat((float) value); }
    public int compareTo(VoxelFloat other) { return value.compareTo(other); }
}
