package edu.jhu.ece.iacl.algorithms.ace;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Anisotropic diffusion calculation for GVF.
 * 
 * @author Blake Lucas
 * 
 */
public class AnisotropicDiffusion extends AbstractCalculation {
	private static final AnisotropicDiffusion ad = new AnisotropicDiffusion();

	public AnisotropicDiffusion() {
		super();
		setLabel("Anisotropic Diffusion");
	}

	public AnisotropicDiffusion(AbstractCalculation parent) {
		super(parent);
		setLabel("Anisotropic Diffusion");
	}

	public static ImageDataFloat doSolve(ImageDataFloat in, ImageDataFloat sphi,
			float alpha, int iterations) {
		return ad.solve(in, sphi, alpha, iterations);
	}

	public ImageDataFloat solve(ImageDataFloat in, ImageDataFloat sphi,
			float alpha, int iterations) {
		return solveVOI(new ImageDataFloat(in), new ImageDataFloat(sphi),
				alpha, iterations);

		/*
		 * int overlap=5+iterations; int XN=in.getRows(); int YN=in.getCols();
		 * int ZN=in.getSlices(); int
		 * sampleRateX=(int)Math.max(1,Math.min(XN/(overlap*3.0),4)); int
		 * sampleRateY=(int)Math.max(1,Math.min(YN/(overlap*3.0),4)); int
		 * sampleRateZ=(int)Math.max(1,Math.min(ZN/(overlap*3.0),4)); int
		 * vxn=XN/sampleRateX; int vyn=YN/sampleRateY; int vzn=ZN/sampleRateZ;
		 * CubicVolumeFloat result=new
		 * CubicVolumeFloat(in.getRows(),in.getCols(),in.getSlices());
		 * CubicVolumeFloat VOIin; CubicVolumeFloat VOIsphi; CubicVolume
		 * VOIresult; int i,j,k,x,y,z; int xv,yv,zv; int beginOffsetX=0; int
		 * endOffsetX=0; int beginOffsetY=0; int endOffsetY=0; int
		 * beginOffsetZ=0; int endOffsetZ=0; int correctOffsetX=0; int
		 * correctOffsetY=0; int correctOffsetZ=0;
		 * 
		 * setTotalUnits((8+iterations)*sampleRateX*sampleRateY*sampleRateZ);
		 * //Caclualte Diffusion for smaller volumes to reduce memory resources
		 * //The volumes are overlapping to allow for correct diffusion
		 * calculation at the boundary for(i=0;i<sampleRateX;i++){
		 * correctOffsetX=i*vxn; beginOffsetX=i*vxn+((i!=0)?-overlap:0);
		 * endOffsetX=(i+1)*vxn+((i!=sampleRateX-1)?overlap:0); for(j=0;j<sampleRateY;j++){
		 * beginOffsetY=j*vyn+((j!=0)?-overlap:0); correctOffsetY=j*vyn;
		 * endOffsetY=(j+1)*vyn+((j!=sampleRateY-1)?overlap:0); for(k=0;k<sampleRateZ;k++){
		 * //setLabel("Anisotropic Diffusion
		 * ("+(i+1)+"/"+sampleRateX+","+(j+1)+"/"+sampleRateY+","+(k+1)+"/"+sampleRateZ+")");
		 * beginOffsetZ=k*vzn+((k!=0)?-overlap:0); correctOffsetZ=k*vzn;
		 * endOffsetZ=(k+1)*vzn+((k!=sampleRateZ-1)?overlap:0); VOIin=new
		 * CubicVolumeFloat(endOffsetX-beginOffsetX,endOffsetY-beginOffsetY,endOffsetZ-beginOffsetZ);
		 * VOIsphi=new
		 * CubicVolumeFloat(endOffsetX-beginOffsetX,endOffsetY-beginOffsetY,endOffsetZ-beginOffsetZ);
		 * for(x=beginOffsetX,xv=0;x<endOffsetX;x++,xv++){
		 * for(y=beginOffsetY,yv=0;y<endOffsetY;y++,yv++){
		 * for(z=beginOffsetZ,zv=0;z<endOffsetZ;z++,zv++){
		 * VOIin.set(xv,yv,zv,inM[x,y,z));
		 * VOIsphi.set(xv,yv,zv,sphiM[x,y,z)); } } }
		 * VOIresult=solveVOI(VOIin,VOIsphi,alpha,iterations); System.gc();
		 * for(x=correctOffsetX,xv=0;xv<vxn;x++,xv++){
		 * for(y=correctOffsetY,yv=0;yv<vyn;y++,yv++){
		 * for(z=correctOffsetZ,zv=0;zv<vzn;z++,zv++){
		 * result.set(x,y,z,VOIresultM[xv+correctOffsetX-beginOffsetX,yv+correctOffsetY-beginOffsetY,zv+correctOffsetZ-beginOffsetZ)); } } } } } }
		 * markCompleted(); return result;
		 */
	}

	private ImageDataFloat solveVOI(ImageDataFloat in,
			ImageDataFloat sphi, double alpha, int iterations) {

		int XN = in.getRows();
		int YN = in.getCols();
		int ZN = in.getSlices();
		int j, k, i, LX, HX, LY, HY, LZ, HZ;
		int index;

		float a11, a12, a13, a22, a23, a33;
		float px, py, pz, mag, px2, py2, pz2;
		float a, b, tmpfloat;
		a = 0.999f;
		b = 0.998f;

		ImageDataFloat S11 = in.mimic();
		ImageDataFloat S12 = in.mimic();
		ImageDataFloat S13 = in.mimic();
		ImageDataFloat S22 = in.mimic();
		ImageDataFloat S23 = in.mimic();
		ImageDataFloat S33 = in.mimic();
		float[][][] S11M = S11.toArray3d();
		float[][][] S12M = S12.toArray3d();
		float[][][] S13M = S13.toArray3d();
		float[][][] S22M = S22.toArray3d();
		float[][][] S23M = S23.toArray3d();
		float[][][] S33M = S33.toArray3d();
		/* Smooth phi *//* Done outside this function */
		float[][][] sphiM = sphi.toArray3d();
		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {

					LX = (i == 0) ? 1 : 0;
					HX = (i == (XN - 1)) ? 1 : 0;
					LY = (j == 0) ? 1 : 0;
					HY = (j == (YN - 1)) ? 1 : 0;
					LZ = (k == 0) ? 1 : 0;
					HZ = (k == (ZN - 1)) ? 1 : 0;
					px = (sphiM[i + 1 - HX][j][k] - sphiM[i - 1 + LX][j][k]) * 0.5f;
					py = (sphiM[i][j + 1 - HY][k] - sphiM[i][j - 1 + LY][k]) * 0.5f;
					pz = (sphiM[i][j][k + 1 - HZ] - sphiM[i][j][k - 1 + LZ]) * 0.5f;

					S11M[i][j][k]=(float)( px * px);
					S12M[i][j][k]=(float)( px * py);
					S13M[i][j][k]=(float)( px * pz);
					S22M[i][j][k]=(float)( py * py);
					S23M[i][j][k]=(float)( py * pz);
					S33M[i][j][k]=(float)( pz * pz);
				}
			}
		}
		incrementCompletedUnits();
		S11 = GaussianBlur(S11, alpha);
		incrementCompletedUnits();
		S12 = GaussianBlur(S12, alpha);
		incrementCompletedUnits();
		S13 = GaussianBlur(S13, alpha);
		incrementCompletedUnits();
		S22 = GaussianBlur(S22, alpha);
		incrementCompletedUnits();
		S23 = GaussianBlur(S23, alpha);
		incrementCompletedUnits();
		S33 = GaussianBlur(S33, alpha);
		incrementCompletedUnits();
		/* Construct the diffusion tensor */
		ImageDataFloat PX2 = in.mimic();
		ImageDataFloat PY2 = in.mimic();
		ImageDataFloat PZ2 = in.mimic();
		float[][][] PX2M = PX2.toArray3d();
		float[][][] PY2M = PY2.toArray3d();
		float[][][] PZ2M = PZ2.toArray3d();

		S11M = S11.toArray3d();
		S12M = S12.toArray3d();
		S13M = S13.toArray3d();
		S22M = S22.toArray3d();
		S23M = S23.toArray3d();
		S33M = S33.toArray3d();
		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {
					a11 = S11M[i][j][k];
					a12 = S12M[i][j][k];
					a13 = S13M[i][j][k];
					a22 = S22M[i][j][k];
					a23 = S23M[i][j][k];
					a33 = S33M[i][j][k];

					/*
					 * Compute the eigenvector corresponging to the largest
					 * eigenvalue, using the Power Method
					 */
					px = 0;
					py = 0;
					pz = 0;
					if (a11 > a22) {
						if (a11 > a33)
							px = 1;
						else
							pz = 1;
					} else {
						if (a22 > a33)
							py = 1;
						else
							pz = 1;
					}

					for (index = 1; index <= 10; index++) {
						px2 = a11 * px + a12 * py + a13 * pz;
						py2 = a12 * px + a22 * py + a23 * pz;
						pz2 = a13 * px + a23 * py + a33 * pz;

						mag = (float) Math.sqrt(px2 * px2 + py2 * py2 + pz2
								* pz2);
						px = px2 / mag;
						py = py2 / mag;
						pz = pz2 / mag;
					}

					/* Now (px, py, pz) is the eigenvector wanted */
					/* Construct the structure tensor */
					S11M[i][j][k] = a - b * px * px;
					S12M[i][j][k] = -b * px * py;
					S13M[i][j][k] = -b * px * pz;
					S22M[i][j][k] = a - b * py * py;
					S23M[i][j][k] = -b * py * pz;
					S33M[i][j][k] = a - b * pz * pz;
				}
			}
		}
		incrementCompletedUnits();
		/* Perform anisotropic diffusion */

		ImageDataFloat out = in.clone();
		float[][][] outM = out.toArray3d();
		for (index = 1; index <= iterations; index++) {
			incrementCompletedUnits();
			for (i = 0; i < XN; i++) {
				for (j = 0; j < YN; j++) {
					for (k = 0; k < ZN; k++) {
						LX = (i == 0) ? 1 : 0;
						LY = (j == 0) ? 1 : 0;
						LZ = (k == 0) ? 1 : 0;

						pz = outM[i][j][k] - outM[i][j][k - 1 + LZ];
						py = outM[i][j][k] - outM[i][j - 1 + LY][k];
						px = outM[i][j][k] - outM[i - 1 + LX][j][k];

						a11 = S11M[i][j][k];
						a12 = S12M[i][j][k];
						a13 = S13M[i][j][k];
						a22 = S22M[i][j][k];
						a23 = S23M[i][j][k];
						a33 = S33M[i][j][k];

						PX2M[i][j][k] = a11 * px + a12 * py + a13 * pz;
						PY2M[i][j][k] = a12 * px + a22 * py + a23 * pz;
						PZ2M[i][j][k] = a13 * px + a23 * py + a33 * pz;
					}
				}
			}

			for (i = 0; i < XN; i++) {
				for (j = 0; j < YN; j++) {
					for (k = 0; k < ZN; k++) {
						HX = (i == (XN - 1)) ? 1 : 0;
						HY = (j == (YN - 1)) ? 1 : 0;
						HZ = (k == (ZN - 1)) ? 1 : 0;

						pz = PZ2M[i][j][k + 1 - HZ] - PZ2M[i][j][k];
						py = PY2M[i][j + 1 - HY][k] - PY2M[i][j][k];
						px = PX2M[i + 1 - HX][j][k] - PX2M[i][j][k];

						tmpfloat = outM[i][j][k] + 0.10f
								* (px + py + pz);
						// if(outM[i,j,k)!=0)if(count++<200)System.out.println(tmpfloat+"
						// "+outM[i,j,k));
						if (tmpfloat < outM[i][j][k])
							tmpfloat = outM[i][j][k];
						if (tmpfloat < 0)
							tmpfloat = 0;
						if (tmpfloat > 1)
							tmpfloat = 1;

						outM[i][j][k] = tmpfloat;
					}
				}
			}
		} /* End of diffusion */
		return out;
	}

	public static ImageDataFloat GaussianBlur(ImageDataFloat input, double alpha) {
		/*
		 * Recursive filtering as in R.Deriche, "Fast Alogorithms for low level
		 * vision", in IEEE PAMI, vol. 12, no. 1, 1990, pp. 78-87. input : the
		 * original image output : smoothed image, can be the same as input
		 * alpha = 1.4/sigma
		 */
		int j, k, i;
		double a1, a2, a3, a4;
		double b0, b1, b2;
		double d1, d2;
		int XN = input.getRows();
		int YN = input.getCols();
		int ZN = input.getSlices();
		ImageDataFloat output = new ImageDataFloat(XN, YN, ZN);
		ImageDataFloat y1 = new ImageDataFloat(XN, YN, ZN);
		ImageDataFloat y2 = new ImageDataFloat(XN, YN, ZN);
		ImageDataFloat r = new ImageDataFloat(XN, YN, ZN);
		float[][][] y1M=y1.toArray3d();
		float[][][] y2M=y2.toArray3d();
		float[][][] rM=r.toArray3d();
		float[][][] inputM=input.toArray3d();
		float[][][] outputM=output.toArray3d();
		d1 = Math.exp(-alpha);
		d2 = d1 * d1;
		b0 = (1.0f - d1) * (1.0f - d1) / (1.0f + 2.0f * alpha * d1 - d2);

		b1 = 2.0f * d1;
		b2 = -d2;

		a1 = b0;
		a2 = b0 * d1 * (alpha - 1.0f);
		a3 = b0 * d1 * (alpha + 1.0f);
		a4 = -b0 * d2;

		/* Smoothing in X-direction */
		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				y1M[i][j][0]=(float)( a1 * inputM[i][j][0]);
				y1M[i][j][1]=(float)( a1 * inputM[i][j][1] + a2
						* inputM[i][j][0] + b1 * y1M[i][j][0]);
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				for (k = 2; k < (int) ZN; k++) {
					y1M[i][j][k]=(float)( a1 * inputM[i][j][k] + a2
							* inputM[i][ j][ k - 1] + b1
							* y1M[i][ j][ k - 1] + b2
							* y1M[i][ j][ k - 2]);
				}
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				y2M[i][j][ZN-1]=(float)( 0);
				y2M[i][j][ZN-2]=(float)( inputM[i][ j][ ZN - 1] * a3 + b1
						* y2M[i][ j][ ZN - 1]);
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				for (k = (int) (ZN - 3); k >= 0; k--) {
					y2M[i][j][k]=(float)( a3 * inputM[i][ j][ k + 1] + a4
							* inputM[i][ j][ k + 2] + b1
							* y2M[i][ j][ k + 1] + b2
							* y2M[i][ j][ k + 2]);
				}
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				for (k = 0; k < (int) ZN; k++) {
					rM[i][j][k]=(float)( y1M[i][j][k] + y2M[i][j][k]);
				}
			}
		}

		/* Smoothing in Y-direction */
		for (i = 0; i < (int) XN; i++) {
			for (k = 0; k < (int) ZN; k++) {
				y1M[i][0][k]=(float)( a1 * rM[i][ 0][ k]);
				y1M[i][1][k]=(float)( a1 * rM[i][ 1][ k] + a2
						* rM[i][ 0][ k] + b1 * y1M[i][ 0][ k]);
			}
		}

		for (i = 0; i < (int) XN; i++) {
			for (k = 0; k < (int) ZN; k++) {
				for (j = 2; j < (int) YN; j++) {
					y1M[i][j][k]=(float)( a1 * rM[i][j][k] + a2
							* rM[i][ j - 1][ k] + b1
							* y1M[i][ j - 1][ k] + b2
							* y1M[i][ j - 2][ k]);
				}
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (k = 0; k < (int) ZN; k++) {
				y2M[i][YN-1][k]=(float)( 0.0);
				y2M[i][YN-2][k]=(float)( a3 * rM[i][ YN - 1][ k] + b1
						* y2M[i][ YN - 1][ k]);
			}
		}
		for (i = 0; i < (int) XN; i++) {
			for (k = 0; k < (int) ZN; k++) {
				for (j = (int) YN - 3; j >= 0; j--) {
					y2M[i][j][k]=(float)( a3 * rM[i][ j + 1][ k] + a4
							* rM[i][ j + 2][ k] + b1
							* y2M[i][ j + 1][ k] + b2
							* y2M[i][ j + 2][ k]);
				}
			}
		}

		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				for (k = 0; k < (int) ZN; k++) {
					rM[i][j][k]=(float)( y1M[i][j][k] + y2M[i][j][k]);
				}
			}
		}
		/* Smoothing in Z-direction */
		for (j = 0; j < (int) YN; j++) {
			for (k = 0; k < (int) ZN; k++) {
				y1M[0][j][k]=(float)( a1 * rM[0][ j][ k]);
				y1M[1][j][k]=(float)( a1 * rM[1][ j][ k] + a2
						* rM[0][ j][ k] + b1 * y1M[0][ j][ k]);
			}
		}

		for (j = 0; j < (int) YN; j++) {
			for (k = 0; k < (int) ZN; k++) {
				for (i = 2; i < (int) XN; i++) {
					y1M[i][j][k]=(float)( a1 * rM[i][j][k] + a2
							* rM[i - 1][ j][ k] + b1
							* y1M[i - 1][ j][ k] + b2
							* y1M[i - 2][ j][ k]);
				}
			}
		}
		for (j = 0; j < (int) YN; j++) {
			for (k = 0; k < (int) ZN; k++) {
				y2M[XN-1][j][k]=(float)( 0.0);
				y2M[XN-2][j][k]=(float)( a3 * rM[XN - 1][ j][ k] + b1
						* y2M[XN - 1][ j][ k]);
			}
		}

		for (j = 0; j < (int) YN; j++) {
			for (k = 0; k < (int) ZN; k++) {
				for (i = (int) XN - 3; i >= 0; i--) {
					y2M[i][j][k]=(float)( a3 * rM[i + 1][ j][ k] + a4
							* rM[i + 2][ j][ k] + b1
							* y2M[i + 1][ j][ k] + b2
							* y2M[i + 2][ j][ k]);
				}
			}
		}

		for (i = 0; i < (int) XN; i++) {
			for (j = 0; j < (int) YN; j++) {
				for (k = 0; k < (int) ZN; k++) {
					outputM[i][j][k]=(float)( y1M[i][j][k]
							+ y2M[i][j][k]);
				}
			}
		}

		return output;
	}
}
