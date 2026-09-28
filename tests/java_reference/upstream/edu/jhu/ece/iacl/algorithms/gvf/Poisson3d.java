package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.structures.image.ImageDataMath;

/**
 * Solve Poisson's equation using multigrid techniques div( g(x) \nabla u(x)) -
 * h(x) (u - f(x)) = 0
 * 
 * @author Blake Lucas
 * 
 */
public class Poisson3d extends AbstractCalculation {
	private static final GVFSolve linSolve = new GVFSolveGaussSeidel();

	private static final DownSample3d downSample = new DownSample3dByHalf();

	private static final Interpolate3d interpolate = new Interpolate3dNearest();

	private static final float MAX_ERROR = 1E-6f;

	private static ImageDataFloat slvsml(ImageDataFloat RHS, ImageDataFloat W, ImageDataFloat HU,
			int maxlevel) {
		ImageDataFloat U = RHS.mimic();
		return linSolve.solve(U, RHS, W, HU, 1, 2, maxlevel);
	}

	private static ImageDataFloat residue(ImageDataFloat U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, int lev, int maxlevel) {
		int i, j, k;
		int h;
		int xn = RHS.getRows();
		int yn = RHS.getCols();
		int zn = RHS.getSlices();
		float cmu;
		ImageDataFloat RES = new ImageDataFloat(xn, yn, zn);
		int prek, nextk, prei, nexti, prej, nextj;
		h = 1 << (maxlevel - lev);
		h = h * h;
		cmu = (float) 1.0 / (float) h;

		float[][][] RESmat=RES.toArray3d();
		float[][][] Umat=U.toArray3d();
		float[][][] RHSmat=RHS.toArray3d();
		float[][][] Wmat=W.toArray3d();
		float[][][] HUmat=HU.toArray3d();
		
		/* If change to half-point symmetric extension, no longer converges */
		for (k = 0; k < xn; k++) {
			prek = (k > 0) ? (k - 1) : 0;
			nextk = (k == (xn - 1)) ? (xn - 1) : (k + 1);
			for (i = 0; i < yn; i++) {
				prei = (i > 0) ? (i - 1) : 0;
				nexti = (i == (yn - 1)) ? (yn - 1) : (i + 1);

				for (j = 0; j < zn; j++) {
					prej = (j > 0) ? (j - 1) : 0;
					nextj = (j == (zn - 1)) ? (zn - 1) : (j + 1);
					RESmat[k][i][j] = HUmat[k][i][j]
							* Umat[k][i][j]
							+ RHSmat[k][i][j]
							- (float)cmu
							* Wmat[k][i][j]
							* (Umat[prek][i][j] + Umat[nextk][i][j]
									+ Umat[k][prei][j] + Umat[k][nexti][j]
									+ Umat[k][i][prej] + Umat[k][i][nextj] - Umat[k][i][j] * 6.0f);
				}
			}
		}
		return RES;
	}

	public Poisson3d() {
		super();
	}

	public Poisson3d(AbstractCalculation parent) {
		super(parent);
	}

	private static Poisson3d poisson3d = new Poisson3d();

	public static ImageDataFloat doSolve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, int ITERS) {
		return poisson3d.solve(_U, RHS, W, HU, ITERS);
	}

	public ImageDataFloat solve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W, ImageDataFloat HU,
			int ITERS) {

		ImageDataFloat U = _U.clone();
		int xn = U.getRows();
		int yn = U.getCols();
		int zn = U.getSlices();

		int m, n, l, mf;
		int iters;
		int n1;
		int maxlevel;

		int j, jcycle, jj;
		int msize, nsize, lsize;

		float maxerr = 10, tmpv;
		ImageDataFloat[] ires;
		ImageDataFloat[] irho;
		ImageDataFloat[] iu;
		ImageDataFloat[] iwei;
		ImageDataFloat[] ih;
		ImageDataFloat[] irhs;
		ImageDataFloat uold = new ImageDataFloat(xn, yn, zn);

		float[][][] Umat=U.toArray3d();
		float[][][] RHSmat=RHS.toArray3d();
		float[][][] Wmat=W.toArray3d();
		float[][][] HUmat=HU.toArray3d();
		
		int[] mo;
		int[] no;
		int[] lo;

		/* Computer maximum level */
		m = 0;
		mf = xn;
		while ((mf >>= 1) > 0)
			m++;
		maxlevel = m;
		m = 0;
		mf = yn;
		while ((mf >>= 1) > 0)
			m++;
		if (maxlevel > m)
			maxlevel = m;
		m = 0;
		mf = zn;
		while ((mf >>= 1) > 0)
			m++;
		if (maxlevel > m)
			maxlevel = m;

		/* Set pre-smoothing steps */
		n1 = 2;

		/* Allocate memory */
		ires = new ImageDataFloat[maxlevel + 1];
		irho = new ImageDataFloat[maxlevel + 1];
		irhs = new ImageDataFloat[maxlevel + 1];
		iu = new ImageDataFloat[maxlevel + 1];
		ih = new ImageDataFloat[maxlevel + 1];
		iwei = new ImageDataFloat[maxlevel + 1];
		mo = new int[maxlevel + 1];
		no = new int[maxlevel + 1];
		lo = new int[maxlevel + 1];

		msize = xn;
		nsize = yn;
		lsize = zn;
		for (j = maxlevel; j >= 1; j--) {
			ires[j] = new ImageDataFloat(msize, nsize, lsize);
			irho[j] = new ImageDataFloat(msize, nsize, lsize);
			irhs[j] = new ImageDataFloat(msize, nsize, lsize);
			iu[j] = new ImageDataFloat(msize, nsize, lsize);
			ih[j] = new ImageDataFloat(msize, nsize, lsize);
			iwei[j] = new ImageDataFloat(msize, nsize, lsize);

			mo[j] = msize;
			no[j] = nsize;
			lo[j] = lsize;
			msize = (msize + 1) >> 1;
			nsize = (nsize + 1) >> 1;
			lsize = (lsize + 1) >> 1;

		}
		for (m = 0; m < xn; m++) {
			for (n = 0; n < yn; n++) {
				for (l = 0; l < zn; l++) {
					iwei[maxlevel].toArray3d()[m][n][l] = Wmat[m][n][l];
					ih[maxlevel].toArray3d()[m][n][l] = HUmat[m][n][l];
					uold.toArray3d()[m][n][l] = 0.0f;
				}
			}
		}

		for (j = maxlevel - 1; j >= 1; j--) {
			iwei[j] = downSample.solve(iwei[j + 1], mo[j + 1], no[j + 1],
					lo[j + 1]);
			ih[j] = downSample
					.solve(ih[j + 1], mo[j + 1], no[j + 1], lo[j + 1]);
		}
		/* Now start multigrid processing */
		setTotalUnits(ITERS);
		for (iters = 1; iters <= ITERS; iters++) {
			irho[maxlevel] = residue(U, RHS, iwei[maxlevel], ih[maxlevel],
					maxlevel, maxlevel);

			maxerr = 0.0f;
			float[][][] irhoArray=irho[maxlevel].toArray3d();
			for (m = 0; m < xn; m++) {
				for (n = 0; n < yn; n++) {
					for (l = 0; l < zn; l++) {
						tmpv = irhoArray[m][n][l];
						if (tmpv < 0)
							tmpv = 0 - tmpv;
						if (maxerr < tmpv)
							maxerr = tmpv;
					}
				}
			}
			if (maxerr < MAX_ERROR)
				break;

			for (j = maxlevel - 1; j >= 1; j--) {

				irho[j] = downSample.solve(irho[j + 1], mo[j + 1], no[j + 1],
						lo[j + 1]);

			}

			iu[1] = slvsml(irho[1], iwei[1], ih[1], maxlevel); /*
																 * Initial
																 * solution on
																 * coarest grid
																 */

			for (j = 2; j <= maxlevel; j++) {

				iu[j] = interpolate.solve(iu[j - 1], iu[j].getRows(), iu[j]
						.getCols(), iu[j].getSlices()); /*
														 * Get the initial guess
														 * of the solution of
														 * the original eq
														 */

				irhs[j] = irho[j].clone();

				for (jcycle = 1; jcycle <= 1; jcycle++) {
					for (jj = j; jj >= 2; jj--) { /* Down stroke of the V */

						iu[jj] = linSolve.solve(iu[jj], irhs[jj], iwei[jj],
								ih[jj], jj, n1, maxlevel);

						ires[jj] = residue(iu[jj], irhs[jj], iwei[jj], ih[jj],
								jj, maxlevel); /* Defect */

						/*
						 * Does the following initialization only need to be
						 * done at jj=2 or doesn't need to be done at all ?? No!
						 * It need to be done at all level, otherwise, the addin
						 * need to be changed! Why not get rid of += in addin,
						 * and get rid of all these initialization?? Am I
						 * missing something?? Indeed, no change should be made
						 * here!
						 */
						iu[jj - 1] = new ImageDataFloat(iu[jj - 1].getRows(),
								iu[jj - 1].getCols(), iu[jj - 1].getSlices());
						irhs[jj - 1] = downSample.solve(ires[jj], mo[jj],
								no[jj], lo[jj]); /*
													 * mf, nf is the size of the
													 * first parameter
													 */

					}

					iu[1] = slvsml(irhs[1], iwei[1], ih[1], maxlevel); /*
																		 * Bottom
																		 * of
																		 * the V
																		 */
					/*
					 * Now iu[1] stores the solution of the error at the coarest
					 * level
					 */

					for (jj = 2; jj <= j; jj++) { /* Upard stroke of V */
						ires[jj] = interpolate.solve(iu[jj - 1], ires[jj]
								.getRows(), ires[jj].getCols(), ires[jj]
								.getSlices());
//						iu[jj].add(ires[jj]);
						ImageDataMath.addFloatImage(iu[jj], ires[jj]);
						/* Post-smoothing */
						iu[jj] = linSolve.solve(iu[jj], irhs[jj], iwei[jj],
								ih[jj], jj, n1, maxlevel);

					}
				}

			}
			maxerr = 0;
//			U.add(iu[maxlevel]);
			ImageDataMath.addFloatImage(U, iu[maxlevel]);
			incrementCompletedUnits();
		} // endof for (iters)
		markCompleted();
		return U;
	}
}
