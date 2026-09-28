package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
/**
 * Generic interface for 3d interpolation
 * @author Blake Lucas
 *
 */
public interface Interpolate3d {
	public ImageDataFloat solve(ImageDataFloat in, int xo, int yo, int zo);
}
