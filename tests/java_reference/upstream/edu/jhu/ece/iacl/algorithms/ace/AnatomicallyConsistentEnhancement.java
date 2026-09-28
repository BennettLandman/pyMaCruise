package edu.jhu.ece.iacl.algorithms.ace;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.algorithms.volume.DistanceFieldAndSkeleton;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.structures.image.ImageDataMath;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.MaskVolume18;

/**
 * Anatomically Consistent Enhancement Algorithm
 * 
 * Creates sulcal boundaries by fast-marching the WM with variable speed. The
 * shock created in the distance function is considered the sulci. The sulci
 * mask is subtracted from the GM to enhance the boundary in sulcal areas.
 * 
 * @author Blake Lucas
 * 
 */
public class AnatomicallyConsistentEnhancement extends AbstractCalculation {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.1 $");
	}


	private ImageData csf;

	private ImageData thinGM;

	private ImageData enhancedGM;
	// private ImageData skel,thinskel,ACEGM;

	private static final int max_iter = 25;

	/**
	 * Solve ACE with intesity threshold = 126.9
	 * 
	 * @param WM
	 *            white matter
	 * @param GM
	 *            gray matter
	 */
	public void solve(ImageData WM, ImageData GM) {
		solve(WM, GM, 126.9);
	}

	/**
	 * get enhanced gray matter
	 * 
	 * @return enhanced gray matter
	 */
	public ImageData getEnhancedGM() {
		return enhancedGM;

	}

	/**
	 * get cerebrospinal fluid
	 * 
	 * @return cerebrospinal fluid
	 */
	public ImageData getSkeleton() {
		return csf;
	}

	/**
	 * get thinned gray matter
	 * 
	 * @return gray matter
	 */
	public ImageData getThinnedSkeleton() {
		return thinGM;
	}

	public AnatomicallyConsistentEnhancement() {
		super();
		setLabel("Anatomically Consistent Enhancement");
	}

	public AnatomicallyConsistentEnhancement(AbstractCalculation calc) {
		super(calc);
		setLabel("Anatomically Consistent Enhancement");
	}

	/**
	 * Solve ACE Resoluts stored in getCSF(), getThinGM() and getEnhancedGM()
	 * 
	 * @param wm
	 *            white matter
	 * @param gm
	 *            gray matter
	 * @param thresh
	 *            intensity threshold
	 */
	public void solve(ImageData WM, ImageData GM, double thresh) {
		ImageDataFloat wm = new ImageDataFloat(WM);
		int rows = wm.getRows();
		int cols = wm.getCols();
		int slices = wm.getSlices();
		int i, j, k, LX, LY, LZ, HX, HY, HZ;
		// int tmpint;
		float tmpint;
		float tmpfloat;
		float px, py, pz, mag;
		// enhancedGM = gm.clone();
		enhancedGM = new ImageDataFloat(GM);
		// ImageData tmp;
		// ImageData skel;
		ImageData closedGM;
		ImageDataFloat PHI = new ImageDataFloat(rows, cols, slices);
		csf = new ImageDataInt(rows, cols, slices);
		ImageDataFloat initcsf = new ImageDataFloat(rows, cols, slices);
		ImageDataFloat tmpVol = new ImageDataFloat(rows, cols, slices);
		setTotalUnits(8);
		DistanceField distField = new DistanceField(this);
		AnisotropicDiffusion diffusion = new AnisotropicDiffusion(this);
		NarrowBandACE narrowBand = new NarrowBandACE(this);
		DistanceFieldAndSkeleton dfas = new DistanceFieldAndSkeleton(this);
		CleanACE clean = new CleanACE(this);
		ThinACE thin = new ThinACE(this);
		// Initialize volume
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// tmpVol.set(i, j, k, thresh - wm.getUByte(i, j, k));
					tmpVol.set(i, j, k, thresh - wm.getFloat(i, j, k));
				}
			}
		}

		PHI = distField.solve(tmpVol, 1000.0f);
		incrementCompletedUnits();
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// if (wm.getUByte(i, j, k) >= 30) {
					if (wm.getFloat(i, j, k) >= 30) {
						tmpVol.set(i, j, k, 1);
					} else {
						tmpVol.set(i, j, k, -1);
					}

					// Borrow tmpvol here

					// tmpint = (int) wm.getUByte(i, j, k)
					// + (int) enhancedGM.getUByte(i, j, k);
					tmpint = wm.getFloat(i, j, k)
							+ enhancedGM.getFloat(i, j, k);
					if (tmpint >= 210) {
						wm.set(i, j, k, 1.0f);
					} else {
						wm.set(i, j, k, 0.0f);
					}

					if (tmpint > 254) {
						csf.set(i, j, k, (int) 0);
					} else {
						csf.set(i, j, k, (int) (254 - (int) tmpint));
					}
				}
			}
		}
		// return type of Morphology is ImageDataInt
		closedGM = Morphology.close(wm, new MaskVolume18(), 4);
		incrementCompletedUnits();
		// int count = 0;
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					initcsf.set(i, j, k, 0);
					// if (closedGM.getUByte(i, j, k) == 1
					// && wm.getUByte(i, j, k) == 0) {
					if (closedGM.getInt(i, j, k) == 1
							&& (int) wm.getFloat(i, j, k) == 0) {
						if (PHI.getFloat(i, j, k) > 0) {
							// if(count++<100)System.out.println("("+i+","+j+","+k+")
							// "+PHI.getFloat(i,j,k)+" "+csf.getByte(i,j,k));
							// tmpfloat = csf.getUByte(i, j, k);
							tmpfloat = csf.getInt(i, j, k);
							initcsf.set(i, j, k, Math.sqrt(tmpfloat / 254.0));
						}
					}
				}
			}
		}
		wm = null;
		// Borrow closegm to record inital phi
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// closedgm = 1 //indicates inside WM isosurface;
					// thus, no ACE is allowed
					if (PHI.getFloat(i, j, k) <= 0
							|| tmpVol.getFloat(i, j, k) > 0)
						closedGM.set(i, j, k, (int) 1);
					else
						closedGM.set(i, j, k, (int) 0);
				}
			}
		}
		tmpVol = AnisotropicDiffusion.GaussianBlur(PHI, 1.8f);

		initcsf = diffusion.solve(initcsf, tmpVol, 1.6f, max_iter);
		incrementCompletedUnits();
		// System.out.println("Begin computing the CSF weighted distance
		// function.\n");
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// Compute the speed function, still store it in initcsf
					initcsf.set(i, j, k, 1.0 - 0.9 * initcsf.getFloat(i, j, k));
				}
			}
		}
		// Blur the speed function
		initcsf = AnisotropicDiffusion.GaussianBlur(initcsf, 1.8f);

		PHI = narrowBand.solve(PHI, (float) 0.1, 1, 10, 5.0);
		incrementCompletedUnits();

		tmpVol = dfas.solve(PHI, initcsf);
		// skel = dfas.getSkeleton();
		incrementCompletedUnits();
		// Compute ACE
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					if (tmpVol.getFloat(i, j, k) <= 0.25 || // reduced from 0.5
							// to 0.25 10-20-03
							// closedGM.getUByte(i, j, k) == 1) {
							closedGM.getInt(i, j, k) == 1) {
						// || phi[k][i][j] > 3.5){ //This condition makes
						// atlas98 failed!
						csf.set(i, j, k, (int) 0);
						continue;
					}

					LX = (i == 0) ? 1 : 0;
					HX = (i == (rows - 1)) ? 1 : 0;
					LY = (j == 0) ? 1 : 0;
					HY = (j == (cols - 1)) ? 1 : 0;
					LZ = (k == 0) ? 1 : 0;
					HZ = (k == (slices - 1)) ? 1 : 0;
					pz = (tmpVol.getFloat(i, j, k + 1 - HZ) - tmpVol.getFloat(
							i, j, k - 1 + LZ)) * 0.5f;
					py = (tmpVol.getFloat(i, j + 1 - HY, k) - tmpVol.getFloat(
							i, j - 1 + LY, k)) * 0.5f;
					px = (tmpVol.getFloat(i + 1 - HX, j, k) - tmpVol.getFloat(i
							- 1 + LX, j, k)) * 0.5f;

					mag = (float) Math.sqrt(px * px + py * py + pz * pz)
							* initcsf.getFloat(i, j, k);
					PHI.set(i, j, k, mag);

					if (mag < 0.85) { // ACE
						tmpfloat = 254.0f * (mag + 0.05f);
						// csf.set(i, j, k, (short) tmpfloat);
						csf.set(i, j, k, (int) tmpfloat);
					} else {
						csf.set(i, j, k, (int) 0); // Indicate non-skeleton
					}
				}
			}
		}
		incrementCompletedUnits();
		clean.solve(csf, 60);
		incrementCompletedUnits();
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// if (csf.getUByte(i, j, k) > 0) {
					if (csf.getInt(i, j, k) > 0) {
						csf.set(i, j, k, (int) 254);
						enhancedGM.set(i, j, k, (float) 0.0f);
					}
				}
			}
		}
		incrementCompletedUnits();
		thinGM = thin.solve(csf);
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					if (i == 0 || j == 0 || k == 0 || i == rows - 1
							|| j == cols - 1 || k == slices - 1) {
						csf.set(i, j, k, (int) 0);
						enhancedGM.set(i, j, k, (float) 0.0f);
					}
				}
			}
		}

		enhancedGM.setHeader(GM.getHeader());
		thinGM.setHeader(GM.getHeader());
		csf.setHeader(GM.getHeader());
		markCompleted();
	}

}
