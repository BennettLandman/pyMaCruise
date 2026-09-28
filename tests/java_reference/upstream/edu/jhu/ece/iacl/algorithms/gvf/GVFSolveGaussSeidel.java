package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Solve system of linear equations using the red/black Gauss Seidel approximation to GVF
 * 
 * mu*del(u)-(u-f_x)*(f_x^2+f_y^2)=0
 * 
 * In the 2d case:
 *u(i,j)=(1-b(i,j)*dt)*u(i,j)+r(u(i+1,j)+u(i,j+1)+u(i-1,j)+u(i,j-1)-4*u(i,j))+c(i,j)*dt
 *
 *b(i,j)=|grad(f(i,j))|
 *c1(i,j)=b(i,j)*f_x(i,j)
 *c2(i,j)=b(i,j)*f_y(i,j)
 *r=mu*dt/(dx*dy)
 *
 * @author Blake Lucas
 * 
 */
public class GVFSolveGaussSeidel implements GVFSolve {
	private static final GVFSolveGaussSeidel linearSolve = new GVFSolveGaussSeidel();

	public static ImageDataFloat doSolve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, int lev, int iter, int maxlevel) {
		return linearSolve.solve(_U, RHS, W, HU, lev, iter, maxlevel);
	}

	public ImageDataFloat solve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W, ImageDataFloat HU,
			int iter, int h) {
		int i;
		float hsquare = h * h;
		ImageDataFloat U = _U.clone();

		for (i = 1; i <= iter; i++) {
			/*
			 * Process black first like follows gives slightly better
			 * convergence
			 */
			U = octant(U, RHS, W, HU, hsquare, 0, 0, 0);
			U = octant(U, RHS, W, HU, hsquare, 0, 1, 1);
			U = octant(U, RHS, W, HU, hsquare, 1, 0, 1);
			U = octant(U, RHS, W, HU, hsquare, 1, 1, 0);
			U = octant(U, RHS, W, HU, hsquare, 1, 0, 0);
			U = octant(U, RHS, W, HU, hsquare, 0, 0, 1);
			U = octant(U, RHS, W, HU, hsquare, 0, 1, 0);
			U = octant(U, RHS, W, HU, hsquare, 1, 1, 1);
		}
		return U;
	}

	public ImageDataFloat solve(ImageDataFloat _U, ImageDataFloat RHS, ImageDataFloat W, ImageDataFloat HU,int lev, int iter, int maxlevel) {
		return solve(_U, RHS, W, HU, iter, 1 << (maxlevel - lev));
	}

	private static ImageDataFloat octant(ImageDataFloat U, ImageDataFloat RHS, ImageDataFloat W,
			ImageDataFloat HU, float hsquare, int kc, int ic, int jc) {

		int x, y, z;
		//CubicVolumeFloat U = _U.clone();
		int xn = U.getRows();
		int yn = U.getCols();
		int zn = U.getSlices();
		int prek, nextk, prei, nexti, prej, nextj;
		float tmpv;
		float[][][] Umat=U.toArray3d();
		float[][][] RHSmat=RHS.toArray3d();
		float[][][] Wmat=W.toArray3d();
		float[][][] HUmat=HU.toArray3d();
		for (x = kc; x < xn; x += 2) {
			prek = (x > 0) ? (x - 1) : 0;
			nextk = (x == (xn - 1)) ? (xn - 1) : (x + 1);
			for (y = ic; y < yn; y += 2) {
				prei = (y > 0) ? (y - 1) : 0;
				nexti = (y == (yn - 1)) ? (yn - 1) : (y + 1);
				for (z = jc; z < zn; z += 2) {
					prej = (z > 0) ? (z - 1) : 0;
					nextj = (z == (zn - 1)) ? (zn - 1) : (z + 1);
					tmpv = HUmat[x][y][z] * hsquare;
					Umat[x][y][z] = ((Umat[prek][y][z] + Umat[nextk][y][z]
							+ Umat[x][prei][z] + Umat[x][nexti][z]
							+ Umat[x][y][prej] + Umat[x][y][nextj])
							* Wmat[x][y][z] - RHSmat[x][y][z] * hsquare)
							/ (6 * Wmat[x][y][z] + tmpv);

				}
			}
		}
		return U;
	}

}
