package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Downsample a volume to dimensions X x Y x Z
 * 
 * @author Blake Lucas
 * 
 */
public abstract class DownSample3d {
	public abstract ImageDataFloat solve(ImageDataFloat in, int x, int y, int z);
}
