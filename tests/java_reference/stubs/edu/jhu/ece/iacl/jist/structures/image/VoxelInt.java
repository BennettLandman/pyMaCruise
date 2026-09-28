package edu.jhu.ece.iacl.jist.structures.image;

/** Minimal integer voxel wrapper for the selected ACE source. */
public final class VoxelInt extends VoxelFloat {
    public VoxelInt() { this(0); }
    public VoxelInt(int value) { super(value); }
}
