package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Nearest Neigbhor interpolation to up sample volume by 2x
 * @author Blake Lucas
 *
 */
public class Interpolate3dNearest implements Interpolate3d {
	private static final Interpolate3dNearest interpolate = new Interpolate3dNearest();

	public static ImageDataFloat doSolve(ImageDataFloat in, int xo, int yo, int zo) {
		return interpolate.solve(in, xo, yo, zo);
	}
	/**
	 * Up sample by 2x
	 * @param in matrix 3d
	 * @param xo ignored
	 * @param yo ignored
	 * @param zo ignored
	 */
	public ImageDataFloat solve(ImageDataFloat in, int xo, int yo, int zo) {
		ImageDataFloat out = new ImageDataFloat(xo, yo, zo);
		int ic, jc, kc, i, j, k;
		// int HX, HY, HZ;
		float[][][] inmat=in.toArray3d();
		float[][][] outmat=out.toArray3d();
		for (kc = 0; kc < xo; kc++) {
			for (ic = 0; ic < yo; ic++) {
				for (jc = 0; jc < zo; jc++) {
					// Down sample by 2
					k = kc >> 1;
					j = jc >> 1;
					i = ic >> 1;
					// HX = (j == (Xi-1))?1:0;
					// HY = (i == (Yi-1))?1:0;
					// HZ = (k==(Zi-1))?1:0;
					// 0.5*(in[k][i][j] + in[k+1-HZ][i+1-HY][j+1-HX]);
					outmat[kc][ic][jc]=inmat[k][i][j];
				}
			}
		}
		return out;
	}
}
