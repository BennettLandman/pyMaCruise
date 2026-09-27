package edu.jhu.ece.iacl.jist.structures.image;

public final class MaskVolume6 {
    public static final int length = 6;
    private static final byte[] X = {1, 0, -1, 0, 0, 0};
    private static final byte[] Y = {0, 1, 0, -1, 0, 0};
    private static final byte[] Z = {0, 0, 0, 0, 1, -1};
    public byte[] getNeighborsX() { return X; }
    public byte[] getNeighborsY() { return Y; }
    public byte[] getNeighborsZ() { return Z; }
}
