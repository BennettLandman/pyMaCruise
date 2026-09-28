package edu.jhu.ece.iacl.algorithms.ace;

import java.util.Vector;

import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelIndexed;
import edu.jhu.ece.iacl.jist.structures.image.VoxelInt;
/**
 * Compute Narrow Band ACE.
 * 
 * @author Blake Lucas
 * 
 */
public class NarrowBandACE extends AbstractCalculation {
	private static final double MaxDist = 6;

	private static final double HowClose = 4;

	private static final NarrowBandACE nb = new NarrowBandACE();

	public NarrowBandACE() {
		super();
		setLabel("Narrow Band ACE");
	}

	public NarrowBandACE(AbstractCalculation parent) {
		super(parent);
		setLabel("Narrow Band ACE");
	}

	public static ImageDataFloat doSolve(ImageDataFloat Phi, double delta_t,
			int maxiters, int L, double KThred) {
		return nb.solve(Phi, delta_t, maxiters, L, KThred);
	}

	public ImageDataFloat solve(ImageDataFloat Phi, double delta_t, int maxiters,
			int L, double KThred) {
		double K; /* Curvature, is the Vs above */
		double Ks = -1.0;
		int j, k, i; /* Index variables */
		double value; /* Temperory storage */
		double tmp; /* Temperory storage */
		int XN = Phi.getRows();
		int YN = Phi.getCols();
		int ZN = Phi.getSlices();
		int lz, ly, lx, hz, hy, hx; /* Varibales implementing mirror boundary */

		/* Variables for implementing Narrow Band */
		int outiter, initer;

		int NBCount; /* Total number of Narrow Band points */
		int LMCount; /* Total number of Landmines */

		int NBIndex, LMIndex; /* Temperory storage */

		int MineHit;

		// FPixel NBPoint;
		// IPixel LMPoint;
		// FPixel *NBStart;
		// IPixel *LMStart;
		// FPixelList NBANDList; /* Narrow Band List */
		// IPixelList LMine; /* LandMine List */

		/*
		 * Dmx (Dmy, Dmz) : backward difference Dpx (Dpy, Dpz) : forward
		 * difference D0x (D0y, D0z) : centered difference
		 */
		double Dmx, Dmy, Dpx, Dpy, D0x, D0y, D0z, Dxx, Dyy, Dxy;
		double Dmz, Dpz, Dzz, Dxz, Dyz;
		double DeltaP; /* Temperory storage */
		double North, South, West, East, Current, Front, Back;
		/* Temperory storage */
		double SD0x, SD0y, SD0z, GPhi; /* Temperory storage */

		KThred = Math.abs(KThred);

		/* Memory allocation */
		ImageDataFloat TempPhi = Phi.mimic();
		VoxelIndexed<VoxelFloat> NBPoint;
		VoxelIndexed<VoxelInt> LMPoint;
		/* Start Narrow Band for Front Propogation */
		Vector<VoxelIndexed<VoxelFloat>> NBANDList = new Vector<VoxelIndexed<VoxelFloat>>();
		Vector<VoxelIndexed<VoxelInt>> LMine = new Vector<VoxelIndexed<VoxelInt>>();
		DistanceField distField=new DistanceField(this);
		setTotalUnits(maxiters);
		for (outiter = 1; outiter <= maxiters; outiter++) {

			/* Find Narrow Band points and Landmines */

			/* Clear old lists */
			NBANDList.clear();
			LMine.clear();
			NBCount = 0;
			LMCount = 0;
			for (i = 0; i < XN; i++) {
				for (j = 0; j < YN; j++) {
					for (k = 0; k < ZN; k++) {
						value = Phi.getDouble(i, j, k);
						if (Math.abs(value) <= MaxDist) {
							/* Put this point to Narrow Band */
							NBPoint = new VoxelIndexed<VoxelFloat>(
									new VoxelFloat((float) value));
							NBPoint.setRefPosition(i, j, k);
							NBANDList.add(NBPoint);
							NBCount++;
							if (Math.abs(value) > HowClose) {
								/* Put it to Landmines */
								LMPoint = new VoxelIndexed<VoxelInt>(
										new VoxelInt((int) Math.signum(value)));
								LMine.add(LMPoint);
							}
						}
					}
				}
			}

			/* Front propogation within Narrow Band */
			MineHit = 0;
			for (initer = 1; initer <= L; initer++) {
				if (MineHit != 0)
					break;
				/* Update Phi at Narrow Band points */
				/*
				 * The new Phi value should first stored at the NBANDList, since
				 * the computation is not parallel
				 */
				for (NBIndex = 1; NBIndex <= NBCount; NBIndex++) {
					NBPoint = NBANDList.get(NBIndex - 1); /*
															 * In case NBCount ==
															 * 0
															 */
					k = NBPoint.getSlice();
					j = NBPoint.getColumn();

					i = NBPoint.getRow();

					ly = (j == 0) ? 1 : 0;
					hy = (j == (YN - 1)) ? 1 : 0;
					lz = (k == 0) ? 1 : 0;
					hz = (k == (ZN - 1)) ? 1 : 0;
					lx = (i == 0) ? 1 : 0;
					hx = (i == (XN - 1)) ? 1 : 0;

					North = Phi.getDouble(i, j - 1 + ly, k);
					South = Phi.getDouble(i, j + 1 - hy, k);
					West = Phi.getDouble(i, j, k - 1 + lz);
					East = Phi.getDouble(i, j, k + 1 - hz);
					Front = Phi.getDouble(i + 1 - hx, j, k);
					Back = Phi.getDouble(i - 1 + lx, j, k);
					Current = Phi.getDouble(i, j, k);

					Dmx = Current - West;
					Dpx = East - Current;
					Dmy = Current - North;
					Dpy = South - Current;
					Dmz = Current - Back;
					Dpz = Front - Current;
					D0x = (East - West) / 2;
					D0y = (South - North) / 2;
					D0z = (Front - Back) / 2;
					Dxx = (West + East - Current - Current);
					Dyy = (North + South - Current - Current);
					Dzz = (Front + Back - Current - Current);
					Dxy = (Phi.getDouble(i, j + 1 - hy, k + 1 - hz)
							+ Phi.getDouble(i, j - 1 + ly, k - 1 + lz)
							- Phi.getDouble(i, j + 1 - hy, k - 1 + lz) - Phi
							.getDouble(i, j - 1 + ly, k + 1 - hz)) / 4;
					Dxz = (Phi.getDouble(i + 1 - hx, j, k + 1 - hz)
							+ Phi.getDouble(i - 1 + lx, j, k - 1 + lz)
							- Phi.getDouble(i + 1 - hx, j, k - 1 + lz) - Phi
							.getDouble(i - 1 + lx, j, k + 1 - hz)) / 4;
					Dyz = (Phi.getDouble(i + 1 - hx, j + 1 - hy, k)
							+ Phi.getDouble(i - 1 + lx, j - 1 + ly, k)
							- Phi.getDouble(i + 1 - hx, j - 1 + ly, k) - Phi
							.getDouble(i - 1 + lx, j + 1 - hy, k)) / 4;
					SD0x = D0x * D0x;
					SD0y = D0y * D0y;
					SD0z = D0z * D0z;
					GPhi = (double) Math.sqrt(SD0x + SD0y + SD0z) + 0.0001;

					K = (Dyy + Dzz)
							* SD0x
							+ (Dxx + Dzz)
							* SD0y
							+ (Dxx + Dyy)
							* SD0z
							- 2
							* (D0x * D0y * Dxy + D0x * D0z * Dxz + D0y * D0z
									* Dyz);
					K = K / (GPhi * GPhi * GPhi);

					/* Threshold K */
					/* if(K > 0) K = 0; */
					if (K > KThred)
						K = KThred;
					else if (K < -KThred)
						K = -KThred;

					/*
					 * K = K/KThred;
					 */

					/* Compute DeltaP */
					tmp = (Dmx > 0) ? Dmx : 0;
					DeltaP = tmp * tmp;
					tmp = (Dpx < 0) ? Dpx : 0;
					DeltaP = DeltaP + tmp * tmp;
					tmp = (Dmy > 0) ? Dmy : 0;
					DeltaP = DeltaP + tmp * tmp;
					tmp = (Dpy < 0) ? Dpy : 0;
					DeltaP = DeltaP + tmp * tmp;
					tmp = (Dmz > 0) ? Dmz : 0;
					DeltaP = DeltaP + tmp * tmp;
					tmp = (Dpz < 0) ? Dpz : 0;
					DeltaP = DeltaP + tmp * tmp;

					DeltaP = (double) Math.sqrt(DeltaP);

					value = DeltaP * 0.01 + 0.95 * Ks * K * GPhi;

					/* Final new value */
					NBANDList.get(NBIndex - 1).set(
							Phi.getDouble(i, j, k) - delta_t * value);
				}

				/* Update Phi */
				for (NBIndex = 1; NBIndex <= NBCount; NBIndex++) {
					NBPoint = (NBANDList.get(NBIndex - 1));
					k = NBPoint.getSlice();
					j = NBPoint.getColumn();
					i = NBPoint.getRow();
					Phi.set(i, j, k, NBPoint.getDouble());
				}

				/* Test LandMines */
				for (LMIndex = 1; LMIndex <= LMCount; LMIndex++) {
					LMPoint = (LMine.get(LMIndex - 1));
					k = LMPoint.getSlice();
					j = LMPoint.getColumn();
					i = LMPoint.getRow();
					if ((int) Math.signum(Phi.getDouble(i, j, k)) != LMPoint
							.getInt()) {
						MineHit = 1;
						break;
					}
				}

			}

			/* Re-initialize using Fast-Marching Method */
			TempPhi = distField.solve(Phi, (float) HowClose);
			Phi = TempPhi.clone();
			incrementCompletedUnits();
		}
		markCompleted();
		return Phi;
	}
}
