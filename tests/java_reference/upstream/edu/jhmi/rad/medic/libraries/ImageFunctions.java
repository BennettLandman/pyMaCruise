package edu.jhmi.rad.medic.libraries;

import java.io.*;
import java.util.*;
import edu.jhmi.rad.medic.utilities.*;

/**
 *
 *  This class computes various basic image functions, e.g.
 *  gradient, convolution, interpolation.
 *	
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ImageFunctions {
	
	// no data: used as a library of functions
	
    // simple image functions

	/**
	 *	minimum value of the image
	 */
    public static float minimum(float[][][] img, int nx, int ny, int nz) {
		float min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static float maximum(float[][][] img, int nx, int ny, int nz) {
		float max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static int minimum(int[][][] img, int nx, int ny, int nz) {
		int min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static int maximum(int[][][] img, int nx, int ny, int nz) {
		int max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static byte minimum(byte[][][] img, int nx, int ny, int nz) {
		byte min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static byte maximum(byte[][][] img, int nx, int ny, int nz) {
		byte max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	normalizes values of the image in [0,1]
	 */
    public static float[] normalize(float[] img, int nx, int ny, int nz, float min, float max) {
		float[] res = new float[nx*ny*nz];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]<min) res[xyz] = min;
			else if (img[xyz]>max) res[xyz] = max;
			else res[xyz] = (img[xyz]-min)/(max-min);
		}
		return res;
	}
	
	/**
	 *	rescale values of the image in [0,1], but allow for values beyond [min,max]
	 */
    public static void rescale(float[] img, int nx, int ny, int nz, float min, float max) {
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			img[xyz] = (img[xyz]-min)/(max-min);
		}
		return;
	}
	
	/**
	 *	rescale values of the image in [0,1], but allow for values beyond [min,max]
	 */
    public static void unscale(float[] img, int nx, int ny, int nz, float min, float max) {
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			img[xyz] = img[xyz]*(max-min) + min;
		}
		return;
	}
	
	/**
	 *	Horn-like derivatives: use a 2x2 forward square as averaging
	 *	the gradient images are passed from above
	 */
    public static void computeGradientHorn(float[][][] img, float[][][] dx, float[][][] dy, float[][][] dz, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=0;x<nx-1;x++) for (y=0;y<ny-1;y++) for (z=0;z<nz-1;z++) {
			dx[x][y][z] = ( (img[x+1][y][z] + img[x+1][y+1][z] + img[x+1][y][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x][y+1][z] + img[x][y][z+1] + img[x][y+1][z+1]) )/4.0f;
			
			dy[x][y][z] = ( (img[x][y+1][z] + img[x+1][y+1][z] + img[x][y+1][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x+1][y][z] + img[x][y][z+1] + img[x+1][y][z+1]) )/4.0f;
			
			dz[x][y][z] = ( (img[x][y][z+1] + img[x+1][y][z+1] + img[x][y+1][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x+1][y][z] + img[x][y+1][z] + img[x+1][y+1][z]) )/4.0f;
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
			dx[x][y][nz-1] = 0.0f;
			dy[x][y][nz-1] = 0.0f;
			dz[x][y][nz-1] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
			dx[x][ny-1][z] = 0.0f;
			dy[x][ny-1][z] = 0.0f;
			dz[x][ny-1][z] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			dx[nx-1][y][z] = 0.0f;
			dy[nx-1][y][z] = 0.0f;
			dz[nx-1][y][z] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple gradient (central differences)
	 */
    public static void computeGradient(float[] img, float[][] J, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			J[0][xyz] = 0.5f*(img[xyz+1] - img[xyz-1]);
			J[1][xyz] = 0.5f*(img[xyz+nx] - img[xyz-nx]);
			J[2][xyz] = 0.5f*(img[xyz+nx*ny] - img[xyz-nx*ny]);
		}
		return;
    }
	
	/**
	 *	simple gradient (central differences)
	 */
    public static void computeGradient(float[][][] img, float[][][] dx, float[][][] dy, float[][][] dz, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=0;x<nx-1;x++) for (y=0;y<ny-1;y++) for (z=0;z<nz-1;z++) {
			dx[x][y][z] = 0.5f*(img[x+1][y][z] - img[x-1][y][z]);
			dy[x][y][z] = 0.5f*(img[x][y+1][z] - img[x][y-1][z]);
			dz[x][y][z] = 0.5f*(img[x][y][z+1] - img[x][y][z-1]);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
			dx[x][y][nz-1] = 0.0f;
			dy[x][y][nz-1] = 0.0f;
			dz[x][y][nz-1] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
			dx[x][ny-1][z] = 0.0f;
			dy[x][ny-1][z] = 0.0f;
			dz[x][ny-1][z] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			dx[nx-1][y][z] = 0.0f;
			dy[nx-1][y][z] = 0.0f;
			dz[nx-1][y][z] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessian(float[][][] img, float[][][][] hessian, int nx, int ny, int nz) {
        int x,y,z,n;
               
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			hessian[x][y][z][0] = img[x+1][y][z] + img[x-1][y][z] - 2.0f*img[x][y][z];
			hessian[x][y][z][1] = img[x][y+1][z] + img[x][y-1][z] - 2.0f*img[x][y][z];
			hessian[x][y][z][2] = img[x][y][z+1] + img[x][y][z-1] - 2.0f*img[x][y][z];
			
			hessian[x][y][z][3] = img[x+1][y+1][z] + img[x-1][y-1][z] - img[x+1][y-1][z] - img[x-1][y+1][z];
			hessian[x][y][z][4] = img[x][y+1][z+1] + img[x][y-1][z-1] - img[x][y+1][z-1] - img[x][y-1][z+1];
			hessian[x][y][z][5] = img[x+1][y][z+1] + img[x-1][y][z-1] - img[x-1][y][z+1] - img[x+1][y][z-1];
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<6;n++) {
			hessian[x][y][0][n] = 0.0f;
			hessian[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<6;n++) {
			hessian[x][0][z][n] = 0.0f;
			hessian[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<6;n++) {
			hessian[0][y][z][n] = 0.0f;
			hessian[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessianEigenvalues(float[][][] img, float[][][][] values, int nx, int ny, int nz) {
        int x,y,z,n;
		float[][] hessian = new float[3][3];
		       
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			hessian[0][0] = img[x+1][y][z] + img[x-1][y][z] - 2.0f*img[x][y][z];
			hessian[1][1] = img[x][y+1][z] + img[x][y-1][z] - 2.0f*img[x][y][z];
			hessian[2][2] = img[x][y][z+1] + img[x][y][z-1] - 2.0f*img[x][y][z];
			
			hessian[0][1] = img[x+1][y+1][z] + img[x-1][y-1][z] - img[x+1][y-1][z] - img[x-1][y+1][z];
			hessian[1][2] = img[x][y+1][z+1] + img[x][y-1][z-1] - img[x][y+1][z-1] - img[x][y-1][z+1];
			hessian[2][0] = img[x+1][y][z+1] + img[x-1][y][z-1] - img[x-1][y][z+1] - img[x+1][y][z-1];
			
			hessian[1][0] = hessian[0][1];
			hessian[2][1] = hessian[1][2];
			hessian[0][2] = hessian[2][0];
			
			// compute eigenvalues
			values[x][y][z] = Numerics.eigenvalues3D(hessian);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<3;n++) {
			values[x][y][0][n] = 0.0f;
			values[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[x][0][z][n] = 0.0f;
			values[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[0][y][z][n] = 0.0f;
			values[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
	 *	gradient shape tensor: return the eigenvalues at each point
	 */
    public static void computeGSTeigenvalues(float[][][] img, float[][][][] values, int nx, int ny, int nz) {
        int x,y,z,n;
		float[][] gst = new float[3][3];
		       
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			// IxIx
			gst[0][0] = 0.25f*(img[x+1][y][z]-img[x-1][y][z])*(img[x+1][y][z]-img[x-1][y][z]);
			// IyIy
			gst[1][1] = 0.25f*(img[x][y+1][z]-img[x][y-1][z])*(img[x][y+1][z]-img[x][y-1][z]);
			// IzIz
			gst[2][2] = 0.25f*(img[x][y][z+1]-img[x][y][z-1])*(img[x][y][z+1]-img[x][y][z-1]);
			// IxIy
			gst[0][1] = 0.25f*(img[x+1][y][z]-img[x-1][y][z])*(img[x][y+1][z]-img[x][y-1][z]);
			gst[1][0] = gst[0][1];
			// IyIz
			gst[1][2] = 0.25f*(img[x][y+1][z]-img[x][y-1][z])*(img[x][y][z+1]-img[x][y][z-1]);
			gst[2][1] = gst[1][2];
			// Izx ?
			gst[2][0] = 0.25f*(img[x][y][z+1]-img[x][y][z-1])*(img[x+1][y][z]-img[x-1][y][z]);
			gst[0][2] = gst[2][0];
			
			// compute eigenvalues
			values[x][y][z] = Numerics.eigenvalues3D(gst);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<3;n++) {
			values[x][y][0][n] = 0.0f;
			values[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[x][0][z][n] = 0.0f;
			values[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[0][y][z][n] = 0.0f;
			values[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
	 *	convolution
	 */
	public static float[][][] convolution(float[][][] image, int nx, int ny, int nz, float[][][] kernel, int kx, int ky, int kz) {
		float[][][] result = new float[nx][ny][nz];
		int xi,yj,zl;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
				xi = x+i; yj = y+j; zl = z+l;
				if ( (xi>=0) && (xi<nx) && (yj>=0) && (yj<ny) && (zl>=0) && (zl<nz) ) {
					result[x][y][z] += image[xi][yj][zl]*kernel[kx+i][ky+j][kz+l];
				}
			}
		}

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	 */
	public static float[][][] separableConvolution(float[][][] image, int nx, int ny, int nz, float[][] kernel, int kx, int ky, int kz) {
		float[][][] result = new float[nx][ny][nz];
		float[][][] temp = new float[nx][ny][nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[x][y][z] += image[x+i][y][z]*kernel[0][kx+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			temp[x][y][z] = 0.0f;	
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[x][y][z] += result[x][y+i][z]*kernel[1][ky+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[x][y][z] += temp[x][y][z+i]*kernel[2][kz+i];
				}
			}
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	 */
	public static float[] separableConvolution(float[] image, int nx, int ny, int nz, float[][] kernel, int kx, int ky, int kz) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[xyz] += image[xyz+i]*kernel[0][kx+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			temp[xyz] = 0.0f;	
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[xyz] += result[xyz+i*nx]*kernel[1][ky+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[xyz] += temp[xyz+i*nx*ny]*kernel[2][kz+i];
				}
			}
		}
		temp = null;

		return result;
	}
		
	/**
	 *	Gaussian kernel
	 */
	public static float[][][] gaussianKernel(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		// create the kernel
		float[][][] kernel = new float[2*kx+1][2*ky+1][2*kz+1];
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
			kernel[kx+i][ky+j][kz+l] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) + (j*j)/(sy*sy) + (l*l)/(sz*sz) );
			sum += kernel[kx+i][ky+j][kz+l];
		}
		// normalize
		for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
			kernel[kx+i][ky+j][kz+l] = kernel[kx+i][ky+j][kz+l]/sum;
		}

		return kernel;
	}
	/**
	 *	Gaussian kernel for separable convolution
	 */
	public static float[][] separableGaussianKernel(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[0] = new float[2*kx+1]; 
		kernel[1] = new float[2*ky+1]; 
		kernel[2] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[0][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[0][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[0][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[1][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[1][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[1][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[2][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[2][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[2][kz+l] /= sum;
		
		return kernel;
	}

	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float nearestNeighborInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return 0.0f;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return 0;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][] image, byte zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][][] image, byte zero, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0][c];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[] image, byte zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static short nearestNeighborInterpolation(short[] image, short zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float[] nearestNeighborInterpolation(float[][] image, float x, float y, float z, int nx, int ny, int nz, int nd) {
		int x0,y0,z0;
		float[] val = new float[3];

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) {
			for (int d=0;d<nd;d++) val[d] = 0.0f;
            return val;
        }
		
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		for (int d=0;d<nd;d++) val[d] = image[d][x0+nx*y0+nx*ny*z0];
		return val;
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static int nearestNeighborInterpolation(int[][][] image, int zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

        return image[x0][y0][z0];
	}
	
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float linearInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[][][] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);
		
		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[][][][] image, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);
		
		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0][c] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0][c] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0][c] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1][c] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0][c] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1][c] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1][c] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1][c];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[][][][] image, float value, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0][c] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0][c] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0][c] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1][c] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0][c] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1][c] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1][c] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1][c];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);

		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float linearInterpolation(float[][][] image, boolean[][][] mask, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val,den;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		den = 0.0f;
		
		if (mask[x0][y0][z0]) {
			val += nalpha*nbeta*ngamma*image[x0][y0][z0];
			den += nalpha*nbeta*ngamma;
		}
		if (mask[x0+1][y0][z0]) {
			val += alpha*nbeta*ngamma*image[x0+1][y0][z0];
			den += alpha*nbeta*ngamma;
		}
		if (mask[x0][y0+1][z0]) {
			val += nalpha*beta*ngamma*image[x0][y0+1][z0];
			den += nalpha*beta*ngamma;
		}
		if (mask[x0][y0][z0+1]) {
			val += nalpha*nbeta*gamma*image[x0][y0][z0+1];
			den += nalpha*nbeta*gamma;
		}
		if (mask[x0+1][y0+1][z0]) {
			val += alpha*beta*ngamma*image[x0+1][y0+1][z0];
			den += alpha*beta*ngamma;
		}
		if (mask[x0][y0+1][z0+1]) {
			val += nalpha*beta*gamma*image[x0][y0+1][z0+1];
			den += nalpha*beta*gamma;
		}
		if (mask[x0+1][y0][z0+1]) {
			val += alpha*nbeta*gamma*image[x0+1][y0][z0+1];
			den += alpha*nbeta*gamma;
		}
		if (mask[x0+1][y0+1][z0+1]) {
			val += alpha*beta*gamma*image[x0+1][y0+1][z0+1];
			den += alpha*beta*gamma;
		}
		if (den>0) {
			return val/den;
		} else {
			return value;
		}
	}
	
	/**
	 *	linear interpolation, with +/- value outside the image
	 *	useful for interpolating gradients (suppose the gradient is +value away from the image)
	 */
	public static float linearGradientInterpolation(float[][][] image, float value, int dim, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with consistent value
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) {
			// dim corresponds to X,Y or Z
				 if ( (dim==0) && (x<0) )    return -value;
			else if ( (dim==0) && (x>nx-2) ) return value;
			else if ( (dim==1) && (y<0) )    return -value;
			else if ( (dim==1) && (y>ny-2) ) return value;
			else if ( (dim==2) && (z<0) )    return -value;
			else if ( (dim==2) && (z>nz-2) ) return value;
			else return 0.0f;
        }
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with distance to boundary outside image
	 */
	public static float linearDistanceInterpolation(float[][][] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;
		float dist = 0.0f;

        // if out of boundary, replace all with consistent value
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) {
			// dim corresponds to X,Y or Z
			if (x<0) dist += (0-x)*value;
			else if (x>nx-2) dist+= (x-nx+2)*value; 
			if (y<0) dist += (0-y)*value;
			else if (y>ny-2) dist+= (y-ny+2)*value; 
			if (z<0) dist += (0-z)*value;
			else if (z>nz-2) dist+= (z-nz+2)*value; 
			return dist;
        }
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 * interpolate a computation mask:
	 * where there is some data, the mask is set to true
	 */
	public static boolean maskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		return inclusiveMaskInterpolation(image, x, y, z, nx, ny, nz);
	}
		
	/**
	 * interpolate a computation mask:
	 * the mask is true if the linear interpolant is above 0.5
	 */
	public static boolean linearMaskInterpolation(boolean[][][] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		if (mask[x0][y0][z0])       val +=nalpha*nbeta*ngamma;
		if (mask[x0+1][y0][z0])     val +=alpha*nbeta*ngamma;
		if (mask[x0][y0+1][z0])     val +=nalpha*beta*ngamma;
		if (mask[x0][y0][z0+1])     val +=nalpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0])   val +=alpha*beta*ngamma;
		if (mask[x0][y0+1][z0+1])   val +=nalpha*beta*gamma;
		if (mask[x0+1][y0][z0+1])   val +=alpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0+1]) val +=alpha*beta*gamma;

		if (val > 0.5) return true;
		else return false;
	}

	/**
	 * interpolate a computation mask:
	 * returns the linearly interpolated value in [0,1]
	 */
	public static float linearMaskInterpolationValue(boolean[][][] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		if (mask[x0][y0][z0])       val +=nalpha*nbeta*ngamma;
		if (mask[x0+1][y0][z0])     val +=alpha*nbeta*ngamma;
		if (mask[x0][y0+1][z0])     val +=nalpha*beta*ngamma;
		if (mask[x0][y0][z0+1])     val +=nalpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0])   val +=alpha*beta*ngamma;
		if (mask[x0][y0+1][z0+1])   val +=nalpha*beta*gamma;
		if (mask[x0+1][y0][z0+1])   val +=alpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0+1]) val +=alpha*beta*gamma;

		return val;
	}

	/**
	 *  interpolate a computation mask:
	 *  where there is some data, the mask is set to true 
	 *  (generates a mask that contain the transformed image)
	 */
	public static boolean inclusiveMaskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        if (image[x0][y0][z0]) return true;
		else if (image[x0+1][y0][z0]) return true; 
        else if (image[x0][y0+1][z0]) return true;
        else if (image[x0][y0][z0+1]) return true; 
        else if (image[x0+1][y0+1][z0]) return true;
        else if (image[x0][y0+1][z0+1]) return true;
        else if (image[x0+1][y0][z0+1]) return true; 
        else if (image[x0+1][y0+1][z0+1]) return true; 

		return false;
	}
	/**
	 *  interpolate a computation mask:
	 *  where there is no data, the mask is set to false
	 *  (generates a mask that is contained in the transformed image)
	 */
	public static boolean exclusiveMaskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        if (!image[x0][y0][z0]) return false;
		else if (!image[x0+1][y0][z0]) return false; 
        else if (!image[x0][y0+1][z0]) return false;
        else if (!image[x0][y0][z0+1]) return false; 
        else if (!image[x0+1][y0+1][z0]) return false;
        else if (!image[x0][y0+1][z0+1]) return false;
        else if (!image[x0+1][y0][z0+1]) return false; 
        else if (!image[x0+1][y0+1][z0+1]) return false; 

		return true;
	}

	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationXderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float beta,gamma,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nbeta*ngamma*(image[x0+1][y0][z0]-image[x0][y0][z0]) 
			+ beta*ngamma*(image[x0+1][y0+1][z0]-image[x0][y0+1][z0])
			+ nbeta*gamma*(image[x0+1][y0][z0+1]-image[x0][y0][z0+1])
			+ beta*gamma*(image[x0+1][y0+1][z0+1]-image[x0][y0+1][z0+1]);

		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationYderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,gamma,nalpha,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*ngamma*(image[x0][y0+1][z0]-image[x0][y0][z0])
			+ alpha*ngamma*(image[x0+1][y0+1][z0]-image[x0+1][y0][z0])
			+ nalpha*gamma*(image[x0][y0+1][z0+1]-image[x0][y0][z0+1])
			+ alpha*gamma*(image[x0+1][y0+1][z0+1]-image[x0+1][y0][z0+1]);
		
		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationZderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,nalpha,nbeta,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		val = nalpha*nbeta*(image[x0][y0][z0+1]-image[x0][y0][z0])
			+ alpha*nbeta*(image[x0+1][y0][z0+1]-image[x0+1][y0][z0])
			+ nalpha*beta*(image[x0][y0+1][z0+1]-image[x0][y0+1][z0])
			+ alpha*beta*(image[x0+1][y0+1][z0+1]-image[x0+1][y0+1][z0]);
		
		return val;
	}
    
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationXderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float beta,gamma,nbeta,ngamma,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		
		xyz0 = x0 + nx*y0 + nx*ny*z0;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nbeta*ngamma*(image[xyz0+1]-image[xyz0]) 
			+ beta*ngamma*(image[xyz0+1+nx]-image[xyz0+nx])
			+ nbeta*gamma*(image[xyz0+1+nx*ny]-image[xyz0+nx*ny])
			+ beta*gamma*(image[xyz0+1+nx+nx*ny]-image[xyz0+nx+nx*ny]);

		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationYderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,gamma,nalpha,ngamma,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		xyz0 = x0 + nx*y0 + nx*ny*z0;

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*ngamma*(image[xyz0+nx]-image[xyz0])
			+ alpha*ngamma*(image[xyz0+1+nx]-image[xyz0+1])
			+ nalpha*gamma*(image[xyz0+nx+nx*ny]-image[xyz0+nx*ny])
			+ alpha*gamma*(image[xyz0+1+nx+nx*ny]-image[xyz0+1+nx*ny]);
		
		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationZderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,nalpha,nbeta,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		xyz0 = x0 + nx*y0 + nx*ny*z0;

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		val = nalpha*nbeta*(image[xyz0+nx*ny]-image[xyz0])
			+ alpha*nbeta*(image[xyz0+1+nx*ny]-image[xyz0+1])
			+ nalpha*beta*(image[xyz0+nx+nx*ny]-image[xyz0+nx])
			+ alpha*beta*(image[xyz0+1+nx+nx*ny]-image[xyz0+1+nx]);
		
		return val;
	}
    
    /**
     *    Setup 3D cubic Lagrangian
     *    @return	the cubic Lagrangian interpolation polynomial kernel
     */
    public static float[][] setup3DCubicLagrangianInterpolation() {
        int i;
        double arg;
		float[][] wt;

        wt = new float[4][1000];
        for ( i = 0; i < 1000; i++ ) {
            arg = i / 999.0;
            wt[0][i] = (float) ( arg * ( 1.0 - arg ) * ( arg - 2.0 ) / 6.0 );
            wt[1][i] = (float) ( ( arg + 1.0 ) * ( arg - 1.0 ) * ( arg - 2.0 ) * 0.5 );
            wt[2][i] = (float) ( arg * ( arg + 1.0 ) * ( 2.0 - arg ) * 0.5 );
            wt[3][i] = (float) ( arg * ( arg + 1.0 ) * ( arg - 1.0 ) / 6.0 );
        }
		return wt;
    }
	
    /**
     *    3D cubic Lagrangian function
	 *	  @param wt		 the interpolation kernel
     *	  @param min	 minimum threshold for interpolated image
     *	  @param max	 maximum threshold for interpolated image
     *    @param x       float point index
     *    @param y       float point index
     *    @param z       float point index
     *    @return        the cubicLagrangian3D interpolated data point
     */
    public static final float cubicLagrangianInterpolation3D(float[][][] image, float [][] wt, float min, float max, float x, float y, float z, int nx, int ny, int nz ) {

        int xbase, ybase, zbase;
        int j0, j1, j2;
        int l0, l1, l2;
        int ix, iy, iz;
        int indexX, indexY, indexZ;
        float diffX, diffY, diffZ;
        float sum;
        float ySum, zSum;
        int offset;

		
		int xD = nx-1;
		int yD = ny-1;
		int zD = nz-1;
		
        xbase = Numerics.floor(x);
        ybase = Numerics.floor(y);
        zbase = Numerics.floor(z);
        diffX = x - xbase;
        diffY = y - ybase;
        diffZ = z - zbase;
        indexX = Numerics.floor( 999.0 * diffX );
        indexY = Numerics.floor( 999.0 * diffY );
        indexZ = Numerics.floor( 999.0 * diffZ );

        // 15% - 20% faster since Math.max and Math.min are function calls
        // I also replaced the Math.abs but saw no speed improvement.
        sum = 0.0f;
        for ( ix = 0, j0 = xbase - 1; j0 <= xbase + 2; ix++, j0++ ) {
            l0 = xD;
            if ( j0 < l0 ) {
                l0 = j0;
            }
            if ( l0 < 0 ) {
                l0 = 0;
            }
            ySum = 0.0f;
            for ( iy = 0, j1 = ybase - 1; j1 <= ybase + 2; iy++, j1++ ) {
                l1 = yD;
                if ( j1 < l1 ) {
                    l1 = j1;
                }
                if ( l1 < 0 ) {
                    l1 = 0;
                }
                zSum = 0.0f;
                for ( iz = 0, j2 = zbase - 1; j2 <= zbase + 2; iz++, j2++ ) {
                    l2 = zD;
                    if ( j2 < l2 ) {
                        l2 = j2;
                    }
                    if ( l2 < 0 ) {
                        l2 = 0;
                    }
                    zSum += wt[iz][indexZ] * image[l0][l1][l2];
                } // for (iz = 0, j2 = zbase - 1; j2 <= zbase + 2;iz++, j2++)
                ySum += wt[iy][indexY] * zSum;
            } // for (iy = 0,j1 = ybase - 1; j1 <= ybase + 2;iy++, j1++)
            sum += wt[ix][indexX] * ySum;
        } // for (ix = 0,j0 = xbase - 1; j0 <= xbase + 2;ix++, j0++)
        if (sum>max) sum = max;
		if (sum<min) sum = min;
		 
        return sum;
    }
    /**
     *    3D cubic Lagrangian function
     *	  @param wt		 the interpolation kernel
     *	  @param min	 minimum threshold for interpolated image
     *	  @param max	 maximum threshold for interpolated image
     *    @param mask    boolean mask: true if the image data is to be used in the interpolation
     *    @param x       float point index
     *    @param y       float point index
     *    @param z       float point index
     *    @return        the cubicLagrangian3D interpolated data point
     */
    public static final float cubicLagrangianInterpolation3D(float[][][] image, boolean[][][] mask, float [][] wt, float min, float max, float x, float y, float z, int nx, int ny, int nz ) {

        int xbase, ybase, zbase;
        int j0, j1, j2;
        int l0, l1, l2;
        int ix, iy, iz;
        int indexX, indexY, indexZ;
        float diffX, diffY, diffZ;
        float sum,den;
        float ySum, zSum, yDen, zDen;
        int offset;

		
		int xD = nx-1;
		int yD = ny-1;
		int zD = nz-1;
		
        xbase = Numerics.floor(x);
        ybase = Numerics.floor(y);
        zbase = Numerics.floor(z);
        diffX = x - xbase;
        diffY = y - ybase;
        diffZ = z - zbase;
        indexX = Numerics.floor( 999.0 * diffX );
        indexY = Numerics.floor( 999.0 * diffY );
        indexZ = Numerics.floor( 999.0 * diffZ );

        // 15% - 20% faster since Math.max and Math.min are function calls
        // I also replaced the Math.abs but saw no speed improvement.
        sum = 0.0f;
		den = 0.0f;
        for ( ix = 0, j0 = xbase - 1; j0 <= xbase + 2; ix++, j0++ ) {
            l0 = xD;
            if ( j0 < l0 ) {
                l0 = j0;
            }
            if ( l0 < 0 ) {
                l0 = 0;
            }
            ySum = 0.0f; yDen = 0.0f;
            for ( iy = 0, j1 = ybase - 1; j1 <= ybase + 2; iy++, j1++ ) {
                l1 = yD;
                if ( j1 < l1 ) {
                    l1 = j1;
                }
                if ( l1 < 0 ) {
                    l1 = 0;
                }
                zSum = 0.0f; zDen = 0.0f;
                for ( iz = 0, j2 = zbase - 1; j2 <= zbase + 2; iz++, j2++ ) {
                    l2 = zD;
                    if ( j2 < l2 ) {
                        l2 = j2;
                    }
                    if ( l2 < 0 ) {
                        l2 = 0;
                    }
					if (mask[l0][l1][l2]) {
						zSum += wt[iz][indexZ] * image[l0][l1][l2];
						zDen += wt[iz][indexZ];
					}
                } // for (iz = 0, j2 = zbase - 1; j2 <= zbase + 2;iz++, j2++)
                ySum += wt[iy][indexY] * zSum;
				yDen += wt[iy][indexY] * zDen;
            } // for (iy = 0,j1 = ybase - 1; j1 <= ybase + 2;iy++, j1++)
            sum += wt[ix][indexX] * ySum;
            den += wt[ix][indexX] * yDen;
        } // for (ix = 0,j0 = xbase - 1; j0 <= xbase + 2;ix++, j0++)
		if (den>0) {
			sum = sum/den;
		} else {
			sum = min;
		}
        if (sum>max) sum = max;
		if (sum<min) sum = min;
		 
        return sum;
    }//cubicLagrangian3D
	
	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[][][] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0][0][0];
		Imax = image[0][0][0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[][][] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0][0][0];
		Imax = image[0][0][0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    GVF-style diffusion on scalar images
     *    @param 	image	the original image, assumed to be in [0,1]
	 *    @param 	ratio	amount of diffusion
	 *    @param 	steps	number of diffusion steps
	 *    @param 	dt		time interval
	 *	  @return 			the diffused image
     */
    public static final float[][][] heatDiffusion(float[][][] image, float ratio, int steps, float dt, int nx, int ny, int nz ) {
		float[][][] diff = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			diff[x][y][z] = image[x][y][z];
		}
		for (int t=0;t<steps;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				diff[x][y][z] += ratio*dt/6.0f*(diff[x-1][y][z]+diff[x+1][y][z]
											   +diff[x][y-1][z]+diff[x][y+1][z]
											   +diff[x][y][z-1]+diff[x][y][z+1]
											   -6.0f*diff[x][y][z])
								 + image[x][y][z]*image[x][y][z]*dt*(image[x][y][z]-diff[x][y][z]);
			}
		}
		return diff;
	}
	/**
     *    GVF-style diffusion on scalar images
     *    @param 	image	the original image
	 *	  @param	scale	maximum image value
	 *    @param 	ratio	amount of diffusion
	 *    @param 	steps	number of diffusion steps
	 *    @param 	dt		time interval
	 *	  @return 			the diffused image
     */
    public static final float[][][] heatDiffusion(float[][][] image, float scale, float ratio, int steps, float dt, int nx, int ny, int nz ) {
		float[][][] diff = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			diff[x][y][z] = image[x][y][z];
		}
		for (int t=0;t<steps;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				diff[x][y][z] += ratio*dt/6.0f*(diff[x-1][y][z]+diff[x+1][y][z]
											   +diff[x][y-1][z]+diff[x][y+1][z]
											   +diff[x][y][z-1]+diff[x][y][z+1]
											   -6.0f*diff[x][y][z])
								 + (image[x][y][z]*image[x][y][z])/(scale*scale)*dt*(image[x][y][z]-diff[x][y][z]);
			}
		}
		return diff;
	}
		
	/**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][][] meanCurvature(float[][][] image, int nx, int ny, int nz ) {
		float[][][] curv = new float[nx][ny][nz];
		float ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		float num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// central differences, no spatial smoothing
			ux = 0.5f*(image[x+1][y][z]-image[x-1][y][z]);
			uy = 0.5f*(image[x][y+1][z]-image[x][y-1][z]);
			uz = 0.5f*(image[x][y][z+1]-image[x][y][z-1]);
				
			uxx = image[x+1][y][z] - 2.0f*image[x][y][z] + image[x-1][y][z];
			uyy = image[x][y+1][z] - 2.0f*image[x][y][z] + image[x][y-1][z];
			uzz = image[x][y][z+1] - 2.0f*image[x][y][z] + image[x][y][z-1];
			
			uxy = 0.25f*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			uyz = 0.25f*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			uzx = 0.25f*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = (float)Math.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0f*den*den*den;
				
			if (den>0) {
				curv[x][y][z] = num/den;
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}
	
	/**
	 *	scale down by a factor
	 */
	public static float[][][] subsample(float[][][] image, int nx, int ny, int nz, int factor) {
		int nsx,nsy,nsz;
		float[][][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx][nsy][nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x][y][z] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x][y][z] += image[x*factor+i][y*factor+j][z*factor+l];
			}
			sub[x][y][z] /= scale;
		}
		return sub;
	}
	/**
	 *	scale down by a factor
	 */
	public static float[] subsample(float[] image, int nx, int ny, int nz, int factor) {
		int nsx,nsy,nsz;
		float[] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x + nsx*y + nsx*nsy*z] += image[x*factor+i + nx*(y*factor+j) + nx*ny*(z*factor+l)];
			}
			sub[x + nsx*y + nsx*nsy*z] /= scale;
		}
		return sub;
	}
	
	/**
	 *	scale down a vector image by a factor
	 */
	public static float[][][][] subsample(float[][][][] image, int nx, int ny, int nz, int nv, int factor) {
		int nsx,nsy,nsz;
		float[][][][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx][nsy][nsz][nv];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int v=0;v<nv;v++) {
			sub[x][y][z][v] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x][y][z][v] += image[x*factor+i][y*factor+j][z*factor+l][v];
			}
			sub[x][y][z][v] /= scale;
		}
		return sub;
	}
	/**
	 *	scale down a vector image by a factor
	 */
	public static float[][] subsample(float[][] image, int nx, int ny, int nz, int nv, int factor) {
		int nsx,nsy,nsz;
		float[][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz][nv];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int v=0;v<nv;v++) {
			sub[x+nsx*y+nsx*nsy*z][v] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x+nsx*y+nsx*nsy*z][v] += image[x*factor+i + nx*(y*factor+j) + nx*ny*(z*factor+l)][v];
			}
			sub[x+nsx*y+nsx*nsy*z][v] /= scale;
		}
		return sub;
	}
	
	
	public static float[][][] surfaceEdgeFilter(float[][][] image, boolean[][][] mask, int nx, int ny, int nz) {
		float[][][] result = new float[nx][ny][nz];
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			result[x][y][z] = 0.0f;
			
			if (mask[x][y][z]) {
				float val, best = 0.0f;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z]-image[x+1][y][z]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z]-image[x+1][y-1][z]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z]-image[x+1][y+1][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y][z-1]-image[x+1][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y][z+1]-image[x+1][y][z+1]
					 +2.0f*image[x][y-1][z-1]-image[x-1][y-1][z-1]-image[x+1][y-1][z-1]
					 +2.0f*image[x][y-1][z+1]-image[x-1][y-1][z+1]-image[x+1][y-1][z+1]
					 +2.0f*image[x][y+1][z-1]-image[x-1][y+1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x][y+1][z+1]-image[x-1][y+1][z+1]-image[x+1][y+1][z+1])/9.0f;			
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z]-image[x][y+1][z]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z]-image[x-1][y+1][z]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z]-image[x+1][y+1][z]
					 +2.0f*image[x][y][z-1]-image[x][y-1][z-1]-image[x][y+1][z-1]
					 +2.0f*image[x][y][z+1]-image[x][y-1][z+1]-image[x][y+1][z+1]
					 +2.0f*image[x-1][y][z-1]-image[x-1][y-1][z-1]-image[x-1][y+1][z-1]
					 +2.0f*image[x-1][y][z+1]-image[x-1][y-1][z+1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y][z-1]-image[x+1][y-1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x+1][y][z+1]-image[x+1][y-1][z+1]-image[x+1][y+1][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y][z-1]-image[x][y][z+1]
					 +2.0f*image[x-1][y][z]-image[x-1][y][z-1]-image[x-1][y][z+1]
					 +2.0f*image[x+1][y][z]-image[x+1][y][z-1]-image[x+1][y][z+1]
					 +2.0f*image[x][y-1][z]-image[x][y-1][z-1]-image[x][y-1][z+1]
					 +2.0f*image[x][y+1][z]-image[x][y+1][z-1]-image[x][y+1][z+1]
					 +2.0f*image[x-1][y-1][z]-image[x-1][y-1][z-1]-image[x-1][y-1][z+1]
					 +2.0f*image[x-1][y+1][z]-image[x-1][y+1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y-1][z]-image[x+1][y-1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x+1][y+1][z]-image[x+1][y+1][z-1]-image[x+1][y+1][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z]-image[x+1][y+1][z]
					 +2.0f*image[x-1][y+1][z]-image[x-2][y][z]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z]-image[x][y-2][z]-image[x+2][y][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y-1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]-image[x][y+2][z-1]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y-1][z+1]-image[x+1][y+1][z+1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]-image[x][y+2][z+1]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z-1]-image[x][y+1][z+1]
					 +2.0f*image[x][y+1][z-1]-image[x][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x][y-1][z+1]-image[x][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]-image[x-1][y+2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]-image[x-1][y][z+2]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]-image[x+1][y+2][z]
					 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]-image[x+1][y][z+2])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z-1]-image[x+1][y][z+1]
					 +2.0f*image[x+1][y][z-1]-image[x][y][z-2]-image[x+2][y][z]
					 +2.0f*image[x-1][y][z+1]-image[x-2][y][z]-image[x][y][z+2]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]-image[x+2][y-1][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]-image[x][y-1][z+2]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]-image[x+2][y+1][z]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]-image[x][y+1][z+2])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y+1][z]-image[x+1][y-1][z]
					 +2.0f*image[x-1][y-1][z]-image[x-2][y][z]-image[x][y-2][z]
					 +2.0f*image[x+1][y+1][z]-image[x][y-2][z]-image[x+2][y][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y+1][z-1]-image[x+1][y-1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]-image[x][y-2][z-1]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y+1][z+1]-image[x+1][y-1][z+1]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]-image[x][y-2][z+1]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z+1]-image[x][y+1][z-1]
					 +2.0f*image[x][y-1][z-1]-image[x][y-2][z]-image[x][y][z-2]
					 +2.0f*image[x][y+1][z+1]-image[x][y][z+2]-image[x][y+2][z]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z+1]-image[x-1][y+1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]-image[x-1][y][z-2]
					 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]-image[x-1][y+2][z]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]-image[x+1][y][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]-image[x+1][y+2][z])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z+1]-image[x+1][y][z-1]
					 +2.0f*image[x-1][y][z-1]-image[x-2][y][z]-image[x][y][z-2]
					 +2.0f*image[x+1][y][z+1]-image[x][y][z+2]-image[x+2][y][z]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z+1]-image[x+1][y-1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]-image[x][y-1][z-2]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]-image[x+2][y-1][z]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]-image[x][y+1][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]-image[x+2][y+1][z])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]-image[x+2][y][z]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]-image[x+2][y+2][z]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]-image[x][y+2][z+2]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]-image[x+2][y][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x+1][y-1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]-image[x-2][y][z]
					 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]-image[x-2][y+2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]-image[x-2][y][z+2]
					 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]-image[x][y+2][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y+1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]-image[x][y][z+2]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]-image[x+2][y][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]-image[x][y-2][z]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]-image[x+2][y-2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]-image[x][y-2][z+2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]-image[x+2][y][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]-image[x+2][y][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]-image[x][y][z-2]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]-image[x+2][y][z-2]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]-image[x][y+2][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]-image[x+2][y+2][z])/7.0f;
				if (val*val>best*best) best = val;
				
				result[x][y][z] = 0.5f*best;
			}
		}
		return result;
	}
	
}
