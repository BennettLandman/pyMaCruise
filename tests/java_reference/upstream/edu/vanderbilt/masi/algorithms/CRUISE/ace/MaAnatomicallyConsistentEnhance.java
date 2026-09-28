package edu.vanderbilt.masi.algorithms.CRUISE.ace;

import java.io.File;
import java.io.IOException;

import  edu.jhu.ece.iacl.algorithms.ace.Morphology;
import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.ace.AnisotropicDiffusion;
import edu.jhu.ece.iacl.algorithms.ace.CleanACE;
import edu.jhu.ece.iacl.algorithms.ace.NarrowBandACE;
import edu.jhu.ece.iacl.algorithms.ace.ThinACE;
import edu.jhu.ece.iacl.algorithms.volume.DistanceFieldAndSkeleton;
import edu.jhu.ece.iacl.jist.io.ImageDataReaderWriter;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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
public class MaAnatomicallyConsistentEnhance extends AbstractCalculation {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.7 $");
	}


	private ImageData csf;

	private ImageData thinGM;

	private ImageData enhancedGM;
	// private ImageData skel,thinskel,ACEGM;

	private static final int max_iter = 25;
	
	// to generate all intermediate files
//	private boolean getIntermediateFiles = true;  

	/**
	 * Solve ACE with intesity threshold = 126.9
	 * 
	 * @param WM
	 *            white matter
	 * @param GM
	 *            gray matter
	 */
	//	public void solve(ImageData WM, ImageData GM) {
	//		solve(WM, GM, 126.9);
	//	}

	/**
	 * get enhanced gray matter
	 * 
	 * @return enhanced gray matter
	 */
	public ImageData getEnhancedGM() {
		return enhancedGM;

	}

	public ImageDataFloat MergeSpeedFunction(ImageDataFloat ImageOrigin, ImageDataFloat AddSpeed, int fMethod){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int rows = ImageOrigin.getRows();
		int cols = ImageOrigin.getCols();
		int slices = ImageOrigin.getSlices();
		int i,j,k;
//		float tempval;
		
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// tmpVol.set(i, j, k, thresh - wm.getUByte(i, j, k));
					if (fMethod == 1){
						if(AddSpeed.getFloat(i, j, k)>0){
							ImageOrigin.set(i, j, k, 1.0f);
						}
					}
					else if(fMethod == 2){
						if(AddSpeed.getFloat(i, j, k)>0){
							ImageOrigin.set(i, j, k, 0.1f);
						}
					}
				}
			}
		}
		
		return ImageOrigin;
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

	public MaAnatomicallyConsistentEnhance() {
		super();
		setLabel("Anatomically Consistent Enhancement");
	}

	public MaAnatomicallyConsistentEnhance(AbstractCalculation calc) {
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
	public void solve(ImageData WM, ImageData GM, double thresh, 
			File dir, ImageData mAcsf, int placeToAddCSF,
			int ReplaceMethod,double alphaG,boolean getIntermediateFiles) {
		ImageDataFloat wm = new ImageDataFloat(WM);
		final ImageDataReaderWriter vrw = ImageDataReaderWriter.getInstance();
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
		ImageHeader refHeader = GM.getHeader();
		// ImageData tmp;
		// ImageData skel;
		ImageData closedGM;
		ImageDataFloat PHI = new ImageDataFloat(rows, cols, slices);
		csf = new ImageDataInt(rows, cols, slices);
		ImageDataFloat initcsf = new ImageDataFloat(rows, cols, slices);
		ImageDataFloat tmpVol = new ImageDataFloat(rows, cols, slices);
		setTotalUnits(8);
		DistanceFieldyh distField = new DistanceFieldyh(this);
		AnisotropicDiffusion diffusion = new AnisotropicDiffusion(this);
		NarrowBandACE narrowBand = new NarrowBandACE(this);
		DistanceFieldAndSkeleton dfas = new DistanceFieldAndSkeleton(this);
		CleanACE clean = new CleanACE(this);
		ThinACE thin = new ThinACE(this);

		if (getIntermediateFiles) {
			try{
				if(!dir.isDirectory()){
					(new File(dir.getCanonicalPath())).mkdir();
				}
			}catch(IOException e){ e.printStackTrace(); }
		}
		// Initialize volume
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					// tmpVol.set(i, j, k, thresh - wm.getUByte(i, j, k));
					tmpVol.set(i, j, k, thresh - wm.getFloat(i, j, k));
				}
			}
		}

		if (getIntermediateFiles) {
			tmpVol.setName("01_tmpVolfirst");
			tmpVol.setHeader(refHeader);
			vrw.write(tmpVol, dir);
		}

		PHI = distField.solve(tmpVol, 1000.0f);

		if (getIntermediateFiles) {
			PHI.setName("02_PHIinit");
			PHI.setHeader(refHeader);
			vrw.write(PHI, dir);
		}

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

		if (getIntermediateFiles) {
			csf.setName("03_csf_naive");
			csf.setHeader(refHeader);
			vrw.write(csf,dir);
		}

		// return type of Morphology is ImageDataInt
		closedGM = Morphology.close(wm, new MaskVolume18(), 4);
		//yk add to show the images
		if (getIntermediateFiles) {
			wm.setName("04_wm_morph");
			vrw.write(wm,dir);
			closedGM.setName("05_closedGM_morph");
			vrw.write(closedGM,dir);
		}

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
		if (getIntermediateFiles) {
			PHI.setName("06_PHI_beforeGaussianBlur");
			vrw.write(PHI,dir);
		}
		tmpVol = AnisotropicDiffusion.GaussianBlur(PHI, 1.8f);
		//yk add to show the images
		if (getIntermediateFiles) {
			tmpVol.setName("07_tmpVol_diffusion");
			tmpVol.setHeader(refHeader);
			vrw.write(tmpVol,dir);
			initcsf.setName("08_initcsf_beforediffusion");
			initcsf.setHeader(refHeader);
			vrw.write(initcsf,dir);
		}

		if (placeToAddCSF==1){
			
			ImageDataFloat initmAcsf =  new ImageDataFloat(mAcsf);
			if (ReplaceMethod==0){
				initcsf = MergeSpeedFunction(initcsf, initmAcsf, 1);
			}else if(ReplaceMethod==1){
				initcsf = initmAcsf;
			}

			
			if (getIntermediateFiles) {
				initcsf.setName("08_5_initcsf_beforediffusion_mAcsf");
				initcsf.setHeader(refHeader);
				vrw.write(initcsf,dir);	
			}
			
		}
		
		initcsf = diffusion.solve(initcsf, tmpVol, 1.6f, max_iter);
		
		
		//		//yk modified to add
		//		File fcsf = new File("/media/huoyuan/2EBA8E41BA8E0615/projects/CRUISE/Test_original_space/TestACE/intermediatefiles/initcsf_beforediffusion_withinterface.nii.gz");
		//		ykinput = vrw.read(fcsf);
		//		ImageDataFloat initcsf_yk =  new ImageDataFloat(ykinput);
		//		initcsf_yk.setName("_07_5_initcsf_beforediffusion_ykmodifed");
		//		vrw.write(initcsf_yk,dir);		
		//		

		if (getIntermediateFiles) {
			initcsf.setName("09_initcsf_diffusion");
			initcsf.setHeader(refHeader);
			vrw.write(initcsf,dir);
		}
		
		if (placeToAddCSF==2){
			ImageDataFloat initmAcsf =  new ImageDataFloat(mAcsf);
			if (ReplaceMethod==0){
				initcsf = MergeSpeedFunction(initcsf, initmAcsf, 1);
			}else if(ReplaceMethod==1){
				initcsf = initmAcsf;
			}
			
			if (getIntermediateFiles) {
				initcsf.setName("09_5_initcsf_afterdiffusion_mAcsf");
				initcsf.setHeader(refHeader);
				vrw.write(initcsf,dir);	
			}

		}
		
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
		if (getIntermediateFiles) {
			initcsf.setName("10_initcsf_speedfunction_beforeGuassianBlur");
			initcsf.setHeader(refHeader);
			vrw.write(initcsf,dir);
		}

		//		initcsf = AnisotropicDiffusion.GaussianBlur(initcsf, 1.8f);
		if (alphaG>0){
			initcsf = AnisotropicDiffusion.GaussianBlur(initcsf, alphaG); //yk modified on 06/13/2015
		}

		if (getIntermediateFiles) {
			initcsf.setName("11_initcsf_GuassianBlur");
			initcsf.setHeader(refHeader);
			vrw.write(initcsf,dir);
			PHI.setName("12_PHI_beforenarrowBand");
			PHI.setHeader(refHeader);
			vrw.write(PHI,dir);
		}
		
		if (placeToAddCSF==3){
			ImageDataFloat initmAcsf =  new ImageDataFloat(mAcsf);
			if (ReplaceMethod==0){
				initcsf = MergeSpeedFunction(initcsf, initmAcsf, 2);
			}else if(ReplaceMethod==1){
				initcsf = initmAcsf;
			}
			if (getIntermediateFiles) {
				initcsf.setName("11_5_initcsf_afterGaussian_mAcsf");
				initcsf.setHeader(refHeader);
				vrw.write(initcsf,dir);	
			}
		}
		

		PHI = narrowBand.solve(PHI, (float) 0.1, 1, 10, 5.0);
		incrementCompletedUnits();

		if (getIntermediateFiles) {
			PHI.setName("13_PHI_narrowBand");
			PHI.setHeader(refHeader);
			vrw.write(PHI,dir);
		}

		tmpVol = dfas.solve(PHI, initcsf);

		if (getIntermediateFiles) {
			tmpVol.setName("14_tmpVol_narrowBand");
			tmpVol.setHeader(refHeader);
			vrw.write(tmpVol,dir);
		}

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

		if (getIntermediateFiles) {
			csf.setName("15_csf_beforeclean");
			vrw.write(csf,dir);
		}

		clean.solve(csf, 60);

		if (getIntermediateFiles) {
			csf.setName("16_csf_afterclean");
			vrw.write(csf,dir);
		}

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

