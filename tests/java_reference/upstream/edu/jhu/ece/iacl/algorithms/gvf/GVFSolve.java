package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
/**
 * Solve system of linear equations
 * @author Blake Lucas
 *
 */
public interface GVFSolve {
	public ImageDataFloat solve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, int lev, int iter, int maxlevel);

	public ImageDataFloat solve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, int iter, int h);
}
