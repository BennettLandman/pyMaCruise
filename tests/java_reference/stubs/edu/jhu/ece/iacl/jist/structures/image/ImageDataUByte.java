package edu.jhu.ece.iacl.jist.structures.image;

public final class ImageDataUByte {
    private final byte[][][] values;
    public ImageDataUByte(int x, int y, int z) { values = new byte[x][y][z]; }
    public byte[][][] toArray3d() { return values; }
    public Byte get(int x, int y, int z) { return values[x][y][z]; }
}
