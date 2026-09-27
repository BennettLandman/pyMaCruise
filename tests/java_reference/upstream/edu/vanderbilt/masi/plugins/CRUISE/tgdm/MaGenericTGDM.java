package edu.vanderbilt.masi.plugins.CRUISE.tgdm;

import java.util.LinkedList;
import java.util.Vector;

import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import edu.jhmi.rad.medic.methods.FastMarching;
import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.graphics.map.SphericalMapCorrection;
import edu.jhu.ece.iacl.algorithms.tgdm.CentralCortexProfile;
import edu.jhu.ece.iacl.algorithms.tgdm.InnerCortexProfile;
import edu.jhu.ece.iacl.algorithms.tgdm.OuterCortexProfile;
import edu.jhu.ece.iacl.algorithms.tgdm.ProfileTGDM;
import edu.jhu.ece.iacl.algorithms.tgdm.ProfileTGDM.Tracking;
import edu.jhu.ece.iacl.algorithms.tgdm.SmoothCortexProfile;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.topology.EighteenConnected;
import edu.jhu.ece.iacl.algorithms.topology.SixConnected;
import edu.jhu.ece.iacl.algorithms.topology.TwentySixConnected;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

public class MaGenericTGDM extends AbstractCalculation {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.5 $");
	}

	// For efficiency reasons, only create one instance of these vectors
	protected Vector3d advectionGrad = new Vector3d();
	protected Vector3d regionGrad = new Vector3d();
	protected Vector3d curvatureGrad = new Vector3d();
	protected Vector3d advectionForce = new Vector3d();

	public static class NBPoint extends Point3i implements Comparable<NBPoint> {
		public double w;// the value
		public int index;

		public int compareTo(NBPoint pt) {
			return (int) Math.signum(Math.abs(this.w) - Math.abs(pt.w));
		}

		public String toString() {
			return "[(" + x + "," + y + "," + z + ") " + w + "]";
		}
	}

	protected NBPoint createNBPoint() {
		return new NBPoint();
	}

	// Test for Simple Point Criterion
	protected static interface SimplePointTest {
		public void test(int d, int i, int j, int[][][] cube, byte[][][] status);
	}

	protected class SimplePointTest18 implements SimplePointTest {
		public void test(int x, int y, int z, int[][][] cube, byte[][][] status) {
			int index;
			int xoff, yoff, zoff;
			int cy, cz, cx;

			for (index = 0; index < 6; index++) {
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cx = x + xoff;
				cy = y + yoff;
				cz = z + zoff;

				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;

				cube[xoff + 1][yoff + 1][zoff + 1] = status[cx][cy][cz];
			}

			for (index = 6; index <= 9; index++) { /* zoff = 0 */
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				if (status[cx][cy][cz] == -1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = -1;
				} else if (cube[0 + 1][yoff + 1][0 + 1] == 1
						|| cube[xoff + 1][0 + 1][0 + 1] == 1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = 1;
				}
			}

			for (index = 10; index <= 13; index++) { /* yoff = 0 */
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				if (status[cx][cy][cz] == -1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = -1;
				} else if (cube[xoff + 1][0 + 1][0 + 1] == 1
						|| cube[0 + 1][0 + 1][zoff + 1] == 1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = 1;
				}
			}

			for (index = 14; index <= 17; index++) { /* xoff = 0 */
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				if (status[cx][cy][cz] == -1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = -1;
				} else if (cube[0 + 1][yoff + 1][0 + 1] == 1
						|| cube[0 + 1][0 + 1][zoff + 1] == 1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = 1;
				}
			}

			for (index = 18; index < 26; index++) {
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				if (status[cx][cy][cz] == -1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = -1;
				} else if (cube[xoff + 1][yoff + 1][0 + 1] == 1
						|| cube[xoff + 1][0 + 1][zoff + 1] == 1
						|| cube[0 + 1][yoff + 1][zoff + 1] == 1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = 1;
				}
			}

		}
	}

	protected class SimplePointTest26 implements SimplePointTest {
		public void test(int x, int y, int z, int[][][] cube, byte[][][] status) {
			/* Test whether this pt is a simple point */
			/* Build geodesic neighborhood, n=26 for inside, 6 outside */
			int index;
			int xoff, yoff, zoff;
			int cy, cz, cx;
			for (index = 0; index < 18; index++) {
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				cube[xoff + 1][yoff + 1][zoff + 1] = status[cx][cy][cz];
			}

			for (index = 18; index < 26; index++) {
				xoff = xoff26[index];
				yoff = yoff26[index];
				zoff = zoff26[index];
				cy = y + yoff;
				cz = z + zoff;
				cx = x + xoff;
				cube[xoff + 1][yoff + 1][zoff + 1] = 0;
				if (cy < 0 || cy >= cols || cz < 0 || cz >= slices || cx < 0
						|| cx >= rows)
					continue;
				if (status[cx][cy][cz] == 1) {
					cube[xoff + 1][yoff + 1][zoff + 1] = 1;
				}
			}
		}
	}

	protected static float DEFAULT_OFFSET = 0.1f;
	protected final static int[] xoff26 = { 0, 1, 0, 0, -1, 0, 1, 1, -1, -1, 1,
			1, -1, -1, 0, 0, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1 };
	protected final static int[] yoff26 = { 1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 0,
			0, 0, 0, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
	protected final static int[] zoff26 = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 1,
			-1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1 };
	protected int rows, cols, slices;
	protected float[][][][] initGrad = null;
	// protected float[][][] initPhi = null;
	protected float[][][] regionForces = null;
	protected float[][][] innerBound = null;
	protected float[][][] outerBound = null;
	protected float[][][] skeletalBound = null;
	protected float[][][] implicitLevelSet = null;
	protected float[][][][] advectionForces = null;
	protected EmbeddedSurface initSurf = null;
	protected EmbeddedSurface finalSurf = null;
	protected ConnectivityRule fgRule;
	protected ConnectivityRule bgRule;
	protected SimplePointTest smpTest;
	protected DistanceField df;
	protected int connectivityRule;
	protected Vector<NBPoint> narrowBand;
	protected LinkedList<NBPoint> landMines;
	protected byte[][][] levelSetSign;
	protected boolean mineHit;
	protected double maxSignedDist;
	protected double timeStep;

	public MaGenericTGDM(AbstractCalculation parent, int rule) {
		super(parent);
		switch (rule) {
		case ConnectivityRule.CONNECT_6_18:
		case ConnectivityRule.CONNECT_18_6:
			this.fgRule = new SixConnected();
			this.bgRule = new EighteenConnected();
			smpTest = new SimplePointTest18();
			break;
		case ConnectivityRule.CONNECT_6_26:
		case ConnectivityRule.CONNECT_26_6:
			this.bgRule = new SixConnected();
			this.fgRule = new TwentySixConnected();
			smpTest = new SimplePointTest26();
			break;
		}
		this.connectivityRule = rule;
		this.setLabel("TGDM");
	}

	public MaGenericTGDM(int rule) {
		super();
		switch (rule) {
		case ConnectivityRule.CONNECT_6_18:
		case ConnectivityRule.CONNECT_18_6:
			this.fgRule = new SixConnected();
			this.bgRule = new EighteenConnected();
			smpTest = new SimplePointTest18();
			break;
		case ConnectivityRule.CONNECT_6_26:
		case ConnectivityRule.CONNECT_26_6:
			this.bgRule = new SixConnected();
			this.fgRule = new TwentySixConnected();
			smpTest = new SimplePointTest26();
			break;
		}
		this.connectivityRule = rule;
		this.setLabel("TGDM");
	}

	protected void endNarrowBand(ProfileTGDM profile) {
		/* Re-initialize using Fast-Marching Method */
		// Use the ObjectProcessing class
		float[][][] TempPhi = df.solve(new ImageDataFloat(implicitLevelSet),
				maxSignedDist).toArray3d();
		/* Make copy */
		int y, x, z; /* Index variables */
		for (x = 0; x < rows; x++) {
			for (y = 0; y < cols; y++) {
				for (z = 0; z < slices; z++) {
					implicitLevelSet[x][y][z] = TempPhi[x][y][z];
					// Outside narrow band
					if ((implicitLevelSet[x][y][z] >= maxSignedDist)
							|| implicitLevelSet[x][y][z] <= (-maxSignedDist))
						continue;
					if (innerBound != null) {
						// Different signs, surfaces do not intersect
						if (implicitLevelSet[x][y][z] * innerBound[x][y][z] < 0)
							continue; /* 0 cannot be ignored */
						// Force level set to same value
						if (implicitLevelSet[x][y][z] > innerBound[x][y][z])
							implicitLevelSet[x][y][z] = innerBound[x][y][z];
					}

					if (outerBound != null) {
						if (implicitLevelSet[x][y][z] * outerBound[x][y][z] < 0)
							continue; /* 0 cannot be ignored */
						if (implicitLevelSet[x][y][z] < outerBound[x][y][z])
							implicitLevelSet[x][y][z] = outerBound[x][y][z];
					}
				}
			}
		}
	}

	protected void endTGDM(ProfileTGDM profile) {
		levelSetSign = null;
		narrowBand.clear();
		landMines.clear();
		IsoSurfaceOnGrid isogrid = new IsoSurfaceOnGrid();
		// Create final surface
		finalSurf = isogrid.solveOriginal(new ImageDataFloat(implicitLevelSet),
				connectivityRule, 0, false);
	}

	protected void execute(ProfileTGDM profile) {
		int maxOuterIters = profile.getOuterIters();
		int maxInnerIters = profile.getInternalIters();
		this.setTotalUnits((maxOuterIters - 1) * profile.getInternalIters()
				+ profile.getFinalInternalIters());
		initTGDM(profile);
		int outerIter, innerIter;
		for (outerIter = 1; outerIter <= maxOuterIters; outerIter++) {
			if (outerIter == maxOuterIters) {
				timeStep = profile.getFinalDeltaT();
				maxInnerIters = profile.getFinalInternalIters();
				maxSignedDist = profile.getFinalMaxSignDist();
			}
			initNarrowBand(profile);
			// Move narrow band
			for (innerIter = 1; innerIter <= maxInnerIters; innerIter++) {
				// landmine hit by surface, must re-initialize narrow band
				if (mineHit)
					break;
				for (NBPoint p : narrowBand) {
					validateUpdatePoint(p, profile, updatePoint(p, profile));
				}
				updateNarrowBand(profile);
				incrementCompletedUnits();
			}
			System.out.println("Completed " + (innerIter - 1)
					+ " iterations out of " + maxInnerIters);
			endNarrowBand(profile);
		}
		endTGDM(profile);
		markCompleted();
	}

	public EmbeddedSurface getFinalSurface() {
		return finalSurf;
	}

	public EmbeddedSurface getInitialSurface() {
		return initSurf;
	}
	
	public ImageDataFloat getFinalLvlSet(){
		return new ImageDataFloat(implicitLevelSet);
	}

	protected void initNarrowBand(ProfileTGDM profile) {
		mineHit = false;
		narrowBand.clear();
		landMines.clear();
		NBPoint nbPoint;
		NBPoint lmPoint;
		double value;
		double MaxDist = profile.getMaxDist();
		double HowClose = profile.getHowClose();
		int index = 0;
		int y, x, z; /* Index variables */
		for (x = 0; x < rows; x++) {
			for (y = 0; y < cols; y++) {
				for (z = 0; z < slices; z++) {
					value = implicitLevelSet[x][y][z];
					// Update sign volume to compute topological numbers
					if (connectivityRule == ConnectivityRule.CONNECT_18_6
							|| connectivityRule == ConnectivityRule.CONNECT_6_26) {
						if (value > 0)
							levelSetSign[x][y][z] = -1; /* Outside */
						else
							levelSetSign[x][y][z] = 1; /* Inside */
					} else {
						if (value > 0)
							levelSetSign[x][y][z] = 1; /* Outside */
						else
							levelSetSign[x][y][z] = -1; /* Inside */
					}
					// distance is less than max distance, create narrow band
					// point
					if (Math.abs(value) <= MaxDist) {
						nbPoint = createNBPoint();
						nbPoint.index = index++;
						nbPoint.x = x;
						nbPoint.y = y;
						nbPoint.z = z;
						nbPoint.w = value;// Use the w field to store the
						// level set value
						narrowBand.add(nbPoint);
						// Point is outside narrow band that needs to be
						// updated, add to land mines
						if (Math.abs(value) > HowClose) {
							/* Put it to Landmines */
							lmPoint = createNBPoint();
							lmPoint.x = x;
							lmPoint.y = y;
							lmPoint.z = z;
							lmPoint.w = SIGN(value);
							landMines.add(lmPoint);
						}
					}
				}
			}

		}
	}

	protected void initTGDM(ProfileTGDM profile) {
		rows = implicitLevelSet.length;
		cols = implicitLevelSet[0].length;
		slices = implicitLevelSet[0][0].length;
		df = new DistanceField(this);
		maxSignedDist = profile.getMaxSignDist();
		// Compute iso-surface from initial level set
		IsoSurfaceOnGrid isogrid = new IsoSurfaceOnGrid();
		initSurf = isogrid.solveOriginal(new ImageDataFloat(implicitLevelSet),
				connectivityRule, 0, false);
		narrowBand = new Vector<NBPoint>();
		landMines = new LinkedList<NBPoint>();

		levelSetSign = new byte[rows][cols][slices];// Notice the different
		// retrieve time step
		// The time step may change during curve evolution
		timeStep = profile.getDeltaT();
	}

	protected int SIGN(double x) {
		if (x > 0)
			return 1;
		else
			return -1;
	}

	public ImageDataFloat solveCentralSurface(ImageData initPhi_,
			ImageDataFloat finalPhi, ImageData VpMImg, ImageData GVFImg_,
			float curvatureForce, float externalForce, float pressureForce,
			int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] initPhi = new float[rows][cols][slices];
		float[][][] VpM = new float[rows][cols][slices];
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					VpM[x][y][z] = (VpMImg.getFloat(x, y, z) - 127.0f) / 127.0f;
					initPhi[x][y][z] = Phi[x][y][z];
				}
			}
		}
		innerBound = initPhi;
		outerBound = (finalPhi != null) ? finalPhi.toArray3d() : null;
		regionForces = VpM;
		implicitLevelSet = Phi;
		skeletalBound = null;
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		execute(new CentralCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		if (advectionForces == null)
			System.err.println("NO Advection Force for Central Surface!");
		System.out.println("External Force " + externalForce);
		System.out.println("Pressure Force " + pressureForce);
		System.out.println("Curvature Force " + curvatureForce);
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(initPhi_.getHeader());
		return retVol;
	}

	public ImageDataFloat solveInnerSurface(ImageData initPhi_,
			ImageDataFloat finalPhi, ImageData VpMImg, ImageData GVFImg_,
			float curvatureForce, float externalForce, float pressureForce,
			float isovalue, float isomem, int maxiters) {
		rows = initPhi_.getRows();
		cols = initPhi_.getCols();
		slices = initPhi_.getSlices();
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] VpM = (VpMImg != null) ? new float[rows][cols][slices]
				: null;
		for (int x = 0; x < rows; x++)
			for (int y = 0; y < cols; y++)
				for (int z = 0; z < slices; z++) {
					Phi[x][y][z] = isovalue - Phi[x][y][z];
				}
		System.out.println("Fast Marching...");
		phi = df.solve(phi, 8);
		Phi = phi.toArray3d();
		// Phi = FastMarching.signedDistanceFunction(Phi, X_Img, Y_Img,
		// Z_Img,8.0f);
		if (VpM != null) {
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
		innerBound = null;
		outerBound = (finalPhi != null) ? finalPhi.toArray3d() : null;
		regionForces = VpM;
		implicitLevelSet = Phi;
		skeletalBound = null;
		this.advectionForces = null;

		execute(new InnerCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		// System.out.println("Finished calcTGDM \n");
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(initPhi_.getHeader());
		return retVol;
	}
	
	
	public ImageDataFloat solveInnerSurfaceFromCentral(ImageData central,
			ImageDataFloat finalPhi, ImageData VpMImg, ImageData GVFImg_,
			float curvatureForce, float externalForce, float pressureForce,
			float isovalue, float isomem, int maxiters) {
		rows = central.getRows();
		cols = central.getCols();
		slices = central.getSlices();
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		
			
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(central);
		float[][][] Phi = phi.toArray3d();
		if (VpMImg != null){
			float[][][] VpM = new float[rows][cols][slices];
			// Phi = FastMarching.signedDistanceFunction(Phi, X_Img, Y_Img,
			// Z_Img,8.0f);
			for (int x = 0; x < rows; x++) {
				for (int y = 0; y < cols; y++) {
					for (int z = 0; z < slices; z++) {
						VpM[x][y][z] = -(VpMImg.getFloat(x, y, z) - isomem * 255)
						/ (isomem * 255.0f);
					}
				}
			}
			regionForces = VpM;
		}else
			regionForces = null ;
		
		System.out.println("Start 3D narrow band \n");
		innerBound = (finalPhi != null) ? finalPhi.toArray3d() : null;
		outerBound = null;
		
		implicitLevelSet = Phi;

		skeletalBound = null;
		

		execute(new InnerCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		// System.out.println("Finished calcTGDM \n");
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(central.getHeader());
		return retVol;
	}

	public ImageDataFloat solveCentralSurfaceFromOuter(ImageData initPhi_,
			ImageData VpMImg, ImageData GVFImg_, float curvatureForce,
			float externalForce, float pressureForce, int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] initPhi = new float[rows][cols][slices];
		float[][][] VpM = new float[rows][cols][slices];
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					VpM[x][y][z] = (VpMImg.getFloat(x, y, z) - 127.0f) / 127.0f;
					initPhi[x][y][z] = Phi[x][y][z];
				}
			}
		}
		innerBound = null;
		outerBound = initPhi;
		regionForces = VpM;
		implicitLevelSet = Phi;
		skeletalBound = null;
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		execute(new CentralCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(initPhi_.getHeader());
		return retVol;
	}

	public ImageDataFloat solveSmoothSurface(ImageData initPhi_,
			float curvatureForce, int maxiters) {
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		rows = initPhi_.getRows();
		cols = initPhi_.getCols();
		slices = initPhi_.getSlices();
		System.out.println("Fast Marching...");
		phi = df.solve(phi, 8);
		Phi = phi.toArray3d();
		float[][][] initPhi = new float[rows][cols][slices];
		double max_dist = 0;
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					float tmp = initPhi[x][y][z] = Phi[x][y][z];
					max_dist = Math.max(max_dist, tmp);
				}
			}
		}
		System.out.println("Initial Max Distance " + max_dist);
		innerBound = null;
		outerBound = null;
		regionForces = null;
		implicitLevelSet = Phi;
		skeletalBound = null;
		this.advectionForces = null;
		regionForces = null;
		execute(new SmoothCortexProfile(-curvatureForce, 0, 0, 0, maxiters));// fphi.toArray3d()
		ImageDataFloat vol = new ImageDataFloat(Phi);
		vol.setName(initPhi_.getName() + "_infvol");
		// initSurf.setName(initPhi_.getName());
		// finalSurf.setName(initPhi_.getName() + "_infsurf");
		vol.setHeader(initPhi_.getHeader());
		return vol;
	}

	public ImageDataFloat solveShrinkWrapSurface(ImageData initPhi_,
			ImageData boundPhi, float curvatureForce, float pressureForce,
			int maxiters) {
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		rows = initPhi_.getRows();
		cols = initPhi_.getCols();
		slices = initPhi_.getSlices();
		System.out.println("Fast Marching...");
		phi = df.solve(phi, 8);
		Phi = phi.toArray3d();
		double max_dist = 0;

		innerBound = (new ImageDataFloat(boundPhi)).toArray3d();
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					float tmp = Phi[x][y][z] = Math.min(innerBound[x][y][z],
							Phi[x][y][z]);
					max_dist = Math.max(max_dist, tmp);
				}
			}
		}
		System.out.println("Initial Max Distance " + max_dist);
		outerBound = null;
		regionForces = null;
		implicitLevelSet = Phi;
		skeletalBound = null;
		this.advectionForces = null;
		regionForces = null;
		execute(new InnerCortexProfile(-curvatureForce, 0, -pressureForce,
				maxiters));// fphi.toArray3d()
		ImageDataFloat vol = new ImageDataFloat(Phi);
		vol.setName(initPhi_.getName() + "_infvol");
		// initSurf.setName(initPhi_.getName());
		// finalSurf.setName(initPhi_.getName() + "_infsurf");
		vol.setHeader(initPhi_.getHeader());
		return vol;
	}

	public ImageDataFloat solveInnerSurfaceFromOuter(ImageData central,
			ImageDataFloat finalPhi, ImageData VpMImg, ImageData GVFImg_,
			float curvatureForce, float externalForce, float pressureForce,
			float isovalue, float isomem, int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(central);
		float[][][] Phi = phi.toArray3d();
		float[][][] VpM = new float[rows][cols][slices];
		// Phi = FastMarching.signedDistanceFunction(Phi, X_Img, Y_Img,
		// Z_Img,8.0f);
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					VpM[x][y][z] = (VpMImg.getFloat(x, y, z) - isomem * 255)
							/ (isomem * 255.0f);
				}
			}
		}
		System.out.println("Start 3D narrow band \n");
		innerBound = null;
		outerBound = (finalPhi != null) ? finalPhi.toArray3d() : null;
		regionForces = VpM;
		implicitLevelSet = Phi;

		skeletalBound = null;
		this.advectionForces = null;

		execute(new InnerCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		// System.out.println("Finished calcTGDM \n");
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(central.getHeader());
		return retVol;
	}

	public ImageDataFloat solveIntermediateSurface(boolean negative,
			ImageData initPhi_, ImageData VpMImg, ImageDataFloat innerVol,
			ImageDataFloat outerVol, float levelSet, float curvatureForce,
			float externalForce, float pressureForce, int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		df = new DistanceField(this);
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] VpM = new float[rows][cols][slices];
		float[][][] img = new float[rows][cols][slices];

		float s;
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					s = VpMImg.getFloat(x, y, z);

					VpM[x][y][z] = ((negative) ? -1 : 1) * (levelSet - s);
				}
			}
		}
		img = FastMarching
				.signedDistanceFunction(VpM, rows, cols, slices, 8.0f);
		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				for (int z = 0; z < slices; z++) {
					VpM[x][y][z] = img[x][y][z] / 8.0f;
				}
			}
		}
		innerBound = (innerVol != null) ? innerVol.toArray3d() : null;
		outerBound = (outerVol != null) ? outerVol.toArray3d() : null;
		regionForces = VpM;
		implicitLevelSet = Phi;
		skeletalBound = null;
		advectionForces = null;
		execute(new SmoothCortexProfile(-curvatureForce, externalForce,
				pressureForce, levelSet, maxiters));// fphi.toArray3d()
		ImageDataFloat vol = new ImageDataFloat(Phi);
		vol.setName("level-" + levelSet);
		vol.setHeader(initPhi_.getHeader());
		return vol;
	}

	public ImageDataFloat solveOuterSurface(ImageData initPhi_,
			ImageDataFloat finalPhi, ImageData VpMImg, ImageData GVFImg_,
			ImageData SkelImg_, float curvatureForce, float externalForce,
			float pressureForce, float isomem, int maxiters) {
		rows = VpMImg.getRows();
		cols = VpMImg.getCols();
		slices = VpMImg.getSlices();
		df = new DistanceField(this);
		float[][][] initPhi = new float[rows][cols][slices];
		ImageDataFloat phi = new ImageDataFloat(initPhi_);
		float[][][] Phi = phi.toArray3d();
		float[][][] VpM = new float[rows][cols][slices];
		ImageDataFloat skel = new ImageDataFloat(SkelImg_);
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
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
		innerBound = initPhi;
		outerBound = (finalPhi != null) ? finalPhi.toArray3d() : null;
		regionForces = VpM;
		implicitLevelSet = Phi;
		skeletalBound = skelvol;
		this.advectionForces = (GVFImg_ != null) ? (new ImageDataFloat(GVFImg_))
				.toArray4d()
				: null;
		execute(new OuterCortexProfile(-curvatureForce, externalForce,
				pressureForce, maxiters));
		ImageDataFloat retVol = new ImageDataFloat(Phi);
		retVol.setHeader(initPhi_.getHeader());
		return retVol;
	}

	protected void updateNarrowBand(ProfileTGDM profile) {
		/* Update Phi */
		int y, x, z; /* Index variables */
		// Copy value from narrow band back to volume, assuming the narrow band
		// value is in the proper range
		for (NBPoint p : narrowBand) {
			x = p.x;
			y = p.y;
			z = p.z;
			implicitLevelSet[x][y][z] = (float) p.w;
			if (innerBound != null) {
				if (p.w * innerBound[x][y][z] < 0)
					continue; /* 0 cannot be ignored */
				if (p.w > innerBound[x][y][z])
					p.w = implicitLevelSet[x][y][z] = innerBound[x][y][z];
			}
			if (outerBound != null) {
				if (p.w * outerBound[x][y][z] < 0)
					continue; /* 0 cannot be ignored */
				if (p.w < outerBound[x][y][z])
					p.w = implicitLevelSet[x][y][z] = outerBound[x][y][z];
			}
		}
		// Test and remember if landmine was hit by any of these update points
		for (NBPoint p : landMines) {
			if (SIGN(implicitLevelSet[p.x][p.y][p.z]) != p.w) {
				mineHit = true;
				break;
			}
		}
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
		double externalWeight = profile.getExternalForce();
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
		// ADVECTION FORCE FUNCTION
		if (advectionForces != null) {
			advectionForce = new Vector3d(advectionForces[x][y][z][0],
					advectionForces[x][y][z][1], advectionForces[x][y][z][2]);
			advectionForce.scale(externalWeight);
			advectionGrad.x = (advectionForce.x > 0) ? Dmx : Dpx;
			advectionGrad.y = (advectionForce.y > 0) ? Dmy : Dpy;
			advectionGrad.z = (advectionForce.z > 0) ? Dmz : Dpz;
		}else{
			advectionForce.scale(0f);
			advectionGrad.scale(0f);
		}

		// System.out.println("UPDATE VALUE "+value);
		double tmp = implicitLevelSet[x][y][z]
				- timeStep
				* profile.getEffectiveForce(K, G, curvatureGrad, regionForce,
						regionGrad, advectionForce, advectionGrad);
		return tmp;
	}

	protected void validateUpdatePoint(NBPoint nbPoint, ProfileTGDM profile,
			double updateVal) {
		int x = nbPoint.x, y = nbPoint.y, z = nbPoint.z;
		int fcon, bcon;
		int[][][] cube = new int[3][3][3];
		if ((implicitLevelSet[x][y][z] > 0 && updateVal <= 0)
				|| (implicitLevelSet[x][y][z] <= 0 && updateVal > 0)) {
			if (skeletalBound != null
					&& (implicitLevelSet[x][y][z] > 0 && skeletalBound[x][y][z] > 0)) {
				nbPoint.w = DEFAULT_OFFSET; /*
											 * Donot allow change
											 */
			} else if (outerBound != null
					&& (implicitLevelSet[x][y][z] > 0 && outerBound[x][y][z] >= 0)) {
				nbPoint.w = DEFAULT_OFFSET; /*
											 * Donot allow change
											 */
			} else if (innerBound != null
					&& (implicitLevelSet[x][y][z] <= 0 && innerBound[x][y][z] <= 0)) {
				nbPoint.w = -DEFAULT_OFFSET; /*
											 * Donot allow change
											 */
			} else {

				smpTest.test(x, y, z, cube, levelSetSign);
				fcon = fgRule.check(cube, 1);
				bcon = bgRule.check(cube, -1);

				if (!(fcon == 1 && bcon == 1)) {

					if (implicitLevelSet[x][y][z] > 0)
						nbPoint.w = DEFAULT_OFFSET;
					else
						nbPoint.w = -DEFAULT_OFFSET;
				} else { // Sign is changed
					nbPoint.w = updateVal;
					levelSetSign[x][y][z] = (byte) -levelSetSign[x][y][z];
				}

			}
		} else { /* No sign change */
			nbPoint.w = updateVal;
		}
	}

}
