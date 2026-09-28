package edu.jhu.ece.iacl.algorithms.gvf;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Down samples a volume by half. The new dimensions are computed automatically
 * 
 * @author Blake Lucas
 * 
 */
public class DownSample3dByHalf extends DownSample3d {
	private static final DownSample3dByHalf downSample = new DownSample3dByHalf();

	public static ImageDataFloat doSolve(ImageDataFloat in) {
		return downSample.solve(in, 0, 0, 0);
	}

	/**
	 * Down sample a volume by half
	 */
	public ImageDataFloat solve(ImageDataFloat in) {
		return solve(in, 0, 0, 0);
	}

	/**
	 * Down sample a volume by half
	 * 
	 * @param in
	 *            3D matrix representation of volume
	 * @param x
	 *            ignored
	 * @param y
	 *            ignored
	 * @param z
	 *            igonred
	 */
	public ImageDataFloat solve(ImageDataFloat in, int x, int y, int z) {
		// System.out.println((in.getRows()<<1)+","+(in.getCols()<<1)+","+(in.getLayers()<<1));
		ImageDataFloat out = new ImageDataFloat((in.getRows() + 1) >> 1,
				(in.getCols() + 1) >> 1, (in.getSlices() + 1) >> 1);

		int xo = out.getRows();
		int yo = out.getCols();
		int zo = out.getSlices();
		int xi = in.getRows();
		int yi = in.getCols();
		int zi = in.getSlices();
		// System.out.println(Zi+" "+Yi+" "+Xi+" "+","+Zo+" "+Yo+" "+Zo);
		int inew, jnew, knew, iold, jold, kold;
		float[][][] outmat=out.toArray3d();
		float[][][] inmat=in.toArray3d();
		for (knew = 0; knew < xo; knew++) {
			for (inew = 0; inew < yo; inew++) {
				for (jnew = 0; jnew < zo; jnew++) {
					kold = knew << 1;
					iold = inew << 1;
					jold = jnew << 1;
					if (kold > (xi - 2))
						kold = xi - 2;
					if (iold > (yi - 2))
						iold = yi - 2;
					if (jold > (zi - 2))
						jold = zi - 2;

					outmat[knew][inew][jnew] = 0.125f * (inmat[kold][iold][jold]
							+ inmat[kold][iold][jold + 1]
							+ inmat[kold][iold + 1][jold]
							+ inmat[kold][iold + 1][jold + 1]
							+ inmat[kold + 1][iold][jold]
							+ inmat[kold + 1][iold][jold + 1]
							+ inmat[kold + 1][iold + 1][jold] + inmat[kold + 1][iold + 1][jold + 1]);

				}
			}
		}
		return out;
	}
}
