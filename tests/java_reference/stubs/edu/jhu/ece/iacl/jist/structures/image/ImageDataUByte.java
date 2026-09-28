package edu.jhu.ece.iacl.jist.structures.image;

public final class ImageDataUByte extends ImageData {
    private final byte[][][] values;
    public ImageDataUByte(int x, int y, int z) {
        values = new byte[x][y][z];
        rows = x; cols = y; slices = z;
    }
    public byte[][][] toArray3d() { return values; }
    public Byte get(int x, int y, int z) { return values[x][y][z]; }
    @Override public int getInt(int x, int y, int z) { return values[x][y][z] & 0xff; }
    @Override public float getFloat(int x, int y, int z) { return values[x][y][z] & 0xff; }
    @Override public double getDouble(int x, int y, int z) { return values[x][y][z] & 0xff; }
    @Override public byte getByte(int x, int y, int z) { return values[x][y][z]; }
    @Override public byte getUByte(int x, int y, int z) { return values[x][y][z]; }
    public void set(int x, int y, int z, byte value) { values[x][y][z] = value; }
}
