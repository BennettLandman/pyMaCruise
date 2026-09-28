package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.gvf.FastMarchingGradient;
import edu.jhu.ece.iacl.algorithms.gvf.FastMarchingGradient.Normalization;
import edu.jhu.ece.iacl.algorithms.tgdm.GenericTGDM.NBPoint;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

public class DynamicPressureTGDM extends GenericTGDM{
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.1 $");
	}


	protected float[][][] curvField=null;
	protected float[][][][] refGradient;
	protected float alpha;
	protected float alignmentForce;
	protected float curvatureBiasForce;
	protected ImageDataFloat referenceLevelSet;
	public float getAlpha() {
		return alpha;
	}
	public void setAlpha(float alpha) {
		this.alpha = alpha;
	}
	public float getAlignmentForce() {
		return alignmentForce;
	}
	public void setAlignmentForce(float alignmentForce) {
		this.alignmentForce = alignmentForce;
	}
	public DynamicPressureTGDM(AbstractCalculation parent, int rule) {
		super(parent, rule);
	}
	public DynamicPressureTGDM(int rule) {
		super(rule);
	}
	public void setReferenceLevelSet(ImageDataFloat levelSet){
		FastMarchingGradient fmg=new FastMarchingGradient();
		innerBound = levelSet.toArray3d();
		this.refGradient=fmg.gradient(innerBound, Normalization.MAGNITUDE);

	}
	public ImageDataFloat solveOuterSurface(ImageData initPhi_,ImageDataFloat finalPhi,
			ImageData VpMImg, ImageData GVFImg_, ImageData SkelImg_,
			float curvatureForce, float externalForce, float pressureForce,
			float isomem, int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		df = new DistanceField(this);
		float[][][] initPhi = new float[rows][cols][slices];
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] VpM = new float[rows][cols][slices];
		ImageDataFloat skel = new ImageDataFloat(SkelImg_);
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(
				GVFImg_)).toArray4d() : null;
				float[][][] skelvol = skel.toArray3d();
				for (int x = 0; x < rows; x++) {
					for (int y = 0; y < cols; y++) {
						for (int z = 0; z < slices; z++) {
							float tmpfloat;
							float vp = VpMImg.getFloat(x, y, z);
							initPhi[x][y][z] = Phi[x][y][z];
							if (Phi[x][y][z] < 0.45f)
								tmpfloat = (vp - 32f) / 32f;

							else if (Phi[x][y][z] > 4.5f)
								tmpfloat = (vp - 220f) / 34f;

							else {
								tmpfloat = (vp - isomem * 254) / (isomem * 254);
								// Boost up force inside CSF
								if (tmpfloat < 0)
									tmpfloat = tmpfloat * 1.5f;
							}

							if (skelvol[x][y][z] > 0)
								tmpfloat = -1.5f;

							if (tmpfloat < -1.5)
								tmpfloat = -1.5f;
							else if (tmpfloat > 1.5)
								tmpfloat = 1.5f;

							VpM[x][y][z] = tmpfloat;
						}
					}
				}
				outerBound = (finalPhi!=null)?finalPhi.toArray3d():null;
				regionForces = VpM;
				implicitLevelSet = Phi;
				skeletalBound = skelvol;
				this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(
						GVFImg_)).toArray4d() : null;
						execute(new OuterCortexProfile(-curvatureForce, externalForce,
								pressureForce, maxiters));
					ImageDataFloat retVol=new ImageDataFloat(Phi);
					retVol.setHeader(initPhi_.getHeader());
					return retVol;
	}
	public ImageDataFloat solveBrainMaskSurface(ImageData initPhi_,ImageDataFloat innerBoundPhi,
			ImageDataFloat VpMImg, ImageData GVFImg_,ImageData imgCurv, float curvatureForce,float externalForce,
			float alignmentForce,float curvatureBiasForce,float alpha,float pressureForce, float isovalue,float isomem,int maxiters) {
		rows = initPhi_.getRows();
		cols = initPhi_.getCols();
		slices = initPhi_.getSlices();
		this.alpha=alpha;
		this.alignmentForce=alignmentForce;
		this.curvatureBiasForce=curvatureBiasForce;
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(
				GVFImg_)).toArray4d() : null;
				ImageDataFloat phi = new ImageDataFloat(initPhi_);
				float[][][] Phi = phi.toArray3d();
				float[][][] VpM =null;
				if(VpMImg!=null){
					VpM=VpMImg.toArray3d();
					for (int x = 0; x < rows; x++) {
						for (int y = 0; y < cols; y++) {
							for (int z = 0; z < slices; z++) {
								VpM[x][y][z] = (VpMImg.getFloat(x, y, z) - isomem * 255)
								/ (isomem * 255.0f);
							}
						}
					}
				}
				System.out.println("Start 3D narrow band \n");
				innerBound=null;
				outerBound = null;
				regionForces = VpM;
				
				implicitLevelSet = Phi;
				skeletalBound = null;
				if(imgCurv!=null){
					System.out.println("Curvature Field "+imgCurv.getName());
					curvField=(new ImageDataFloat(imgCurv)).toArray3d();
				}
				this.advectionForces = null;
				DynamicPressureCortexProfile profile=new DynamicPressureCortexProfile(-curvatureForce, externalForce,
						pressureForce, maxiters);
				execute(profile);
				// System.out.println("Finished calcTGDM \n");
				ImageDataFloat retVol=new ImageDataFloat(Phi);
				retVol.setHeader(initPhi_.getHeader());
				return retVol;
	}
	protected double updatePoint(NBPoint nbPoint, ProfileTGDM profile) {
		int LX, LY, LZ, HX, HY, HZ; /* Varibales implementing mirror boundary */
		double North, South, West, East, Current, Front, Back;
		/*
		 * Dmx (Dmy, Dmz) : backward difference Dpx (Dpy, Dpz) : forward
		 * difference D0x (D0y, D0z) : centered difference
		 */
		double Dmx, Dmy, Dpx, Dpy, D0x, D0y, D0z, Dxx, Dyy, Dxy;
		double Dmz, Dpz, Dzz, Dxz, Dyz;
		double SD0x, SD0y, SD0z; /* Temperory storage */
		int x, y, z;
		double K, G; /* Curvature, is the Vs above */

		double regionWeight = profile.getPressureForce();
		double regionForce; /* pressure force */
		x = nbPoint.x;
		y = nbPoint.y;
		z = nbPoint.z;
		LY = (y == 0) ? 1 : 0;
		HY = (y == (cols - 1)) ? 1 : 0;
		LZ = (z == 0) ? 1 : 0;
		HZ = (z == (slices - 1)) ? 1 : 0;
		LX = (x == 0) ? 1 : 0;
		HX = (x == (rows - 1)) ? 1 : 0; 
		
		North = implicitLevelSet[x][y - 1 + LY][z];
		South = implicitLevelSet[x][y + 1 - HY][z];
		West = implicitLevelSet[x - 1 + LX][y][z];
		East = implicitLevelSet[x + 1 - HX][y][z];
		Front = implicitLevelSet[x][y][z + 1 - HZ];
		Back = implicitLevelSet[x][y][z - 1 + LZ];
		Current = implicitLevelSet[x][y][z];
		// m for minus
		// p for plus
		Dmx = Current - West;
		Dpx = East - Current;
		Dmy = Current - North;
		Dpy = South - Current;
		Dmz = Current - Back;
		Dpz = Front - Current;
		// central differences
		D0x = (East - West) / 2;
		D0y = (South - North) / 2;
		D0z = (Front - Back) / 2;
		// second derivative
		Dxx = (West + East - Current - Current);
		Dyy = (North + South - Current - Current);
		Dzz = (Front + Back - Current - Current);
		// cross derivatives
		Dxy = (implicitLevelSet[x + 1 - HX][y + 1 - HY][z]
		                                                + implicitLevelSet[x - 1 + LX][y - 1 + LY][z]
		                                                                                           - implicitLevelSet[x - 1 + LX][y + 1 - HY][z] - implicitLevelSet[x
		                                                                                                                                                            + 1 - HX][y - 1 + LY][z]) / 4;

		Dxz = (implicitLevelSet[x + 1 - HX][y][z + 1 - HZ]
		                                       + implicitLevelSet[x - 1 + LX][y][z - 1 + LZ]
		                                                                         - implicitLevelSet[x - 1 + LX][y][z + 1 - HZ] - implicitLevelSet[x
		                                                                                                                                          + 1 - HX][y][z - 1 + LZ]) / 4;

		Dyz = (implicitLevelSet[x][y + 1 - HY][z + 1 - HZ]
		                                       + implicitLevelSet[x][y - 1 + LY][z - 1 + LZ]
		                                                                         - implicitLevelSet[x][y - 1 + LY][z + 1 - HZ] - implicitLevelSet[x][y
		                                                                                                                                             + 1 - HY][z - 1 + LZ]) / 4;
		// PHI GRADIENT COMPONENTS SQUARED
		SD0x = D0x * D0x;
		SD0y = D0y * D0y;
		SD0z = D0z * D0z;

		// MEAN CURVATURE NUMERATOR
		K = (Dyy + Dzz) * SD0x + (Dxx + Dzz) * SD0y + (Dxx + Dyy) * SD0z - 2
		* (D0x * D0y * Dxy + D0x * D0z * Dxz + D0y * D0z * Dyz);
		// GAUSSIAN CURVATURE NUMERATOR
		G = (Dyy * Dzz - Dyz * Dyz)
		* SD0x
		+ (Dxx * Dzz - Dxz * Dxz)
		* SD0y
		+ (Dxx * Dyy - Dxy * Dxy)
		* SD0z
		+ 2
		* (D0x * D0y * (Dxz * Dyz - Dxy * Dzz) + D0x * D0z
				* (Dxy * Dyz - Dxz * Dyy) + D0y * D0z
				* (Dxy * Dxz - Dyz * Dxx));
		// CENTRAL DIFFERENCE PHI GRADIENT MAGNITUDE

		// Variable pressure field
		curvatureGrad.x = D0x;
		curvatureGrad.y = D0y;
		curvatureGrad.z = D0z;
		if (regionForces == null) {
			regionForce = regionWeight;
		} else {
			regionForce = regionWeight * regionForces[x][y][z];
		}
		
		if (regionForce > 0) {
			regionGrad.x = ((Dmx > 0) ? Dmx : 0) + ((Dpx < 0) ? Dpx : 0);
			regionGrad.y = ((Dmy > 0) ? Dmy : 0) + ((Dpy < 0) ? Dpy : 0);
			regionGrad.z = ((Dmz > 0) ? Dmz : 0) + ((Dpz < 0) ? Dpz : 0);
		} else {
			regionGrad.x = ((Dmx < 0) ? Dmx : 0) + ((Dpx > 0) ? Dpx : 0);
			regionGrad.y = ((Dmy < 0) ? Dmy : 0) + ((Dpy > 0) ? Dpy : 0);
			regionGrad.z = ((Dmz < 0) ? Dmz : 0) + ((Dpz > 0) ? Dpz : 0);
		}
		if(refGradient!=null){
			double len=regionGrad.length();
			//Add force to regularize thickness my minimizing inner product between inner and outer level set gradient normals. 
			double angularForce=0.5-0.5*(refGradient[x][y][z][0]*regionGrad.x+refGradient[x][y][z][1]*regionGrad.y+refGradient[x][y][z][2]*regionGrad.z)/((len>0)?len:1);
			if(alpha!=1){
				angularForce=Math.pow(Math.abs(angularForce),alpha);
			}
			//Add pressure force to push surface outwards in regions that are biased towards higher curvature.
			if(curvField!=null){
				regionForce=regionForce-alignmentForce*angularForce+curvField[x][y][z]*curvatureBiasForce;
			} else {
				regionForce=regionForce-alignmentForce*angularForce;				
			}
		}

		// System.out.println("UPDATE VALUE "+value);
		double tmp = implicitLevelSet[x][y][z]
		                                    - timeStep
		                                    * profile.getEffectiveForce(K, G, curvatureGrad, regionForce,
		                                    		regionGrad, advectionForce, advectionGrad);
		return tmp;
	}

}
