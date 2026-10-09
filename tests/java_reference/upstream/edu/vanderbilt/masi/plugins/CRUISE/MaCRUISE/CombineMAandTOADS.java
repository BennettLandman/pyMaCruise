package edu.vanderbilt.masi.plugins.CRUISE.MaCRUISE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFileCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.vanderbilt.masi.plugins.CRUISE.utilities.MASImorph;

/**
 * @author Yuankai Huo
 * @email yuankai.huo@vanderbilt.edu
 * @version 1.0
 *
 *transfer the toads cruise space image to original space
 *
 */

public class CombineMAandTOADS extends ProcessingAlgorithm{

	// multi-atlas segmentations
	private ParamVolume GM_MA,WM_MA,BrainMask_MA;
	// TOADS segmentations
	private ParamVolume GM_toads,WM_toads,Mask_toads;
	// combined segmentations
	private ParamVolume GM_combine,WM_combine,Mask_combine;
	
	private ParamOption ErodeMethod;
	private ParamDouble GM_confident_level; 

	private static final String cvsversion = "$Revision: 1.2 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "combine multi-atlas segmentations and TOADS segmentations" + revnum + "\n";

	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.setPackage("Vanderbilt");
		inputParams.setCategory("format");
		inputParams.setLabel("Combine MA and TOADS");
		inputParams.setName("Combine_MA_and_TOADS");

		// input multi-atlas segmentation
		GM_MA=new ParamVolume("multiatlas_GM volume");
		WM_MA=new ParamVolume("multiatlas_WM volume");
		BrainMask_MA=new ParamVolume("multiatlas_BrainMask volume");
		
		// TOADS segmentations
		GM_toads=new ParamVolume("toads_GM volume");
		WM_toads=new ParamVolume("toads_WM volume");
		Mask_toads=new ParamVolume("toads_Mask volume");
		
		// 3D min filter or slice by slice like matlab
		ErodeMethod = new ParamOption("Erode Method", new String[] {
				"3D", "3Dmatlab"});
		ErodeMethod.setValue(0);
		ErodeMethod.setMandatory(false);
		
		// define the difference between 
		GM_confident_level = new ParamDouble("GMconfident level");
		GM_confident_level.setValue(0.001);
		GM_confident_level.setMandatory(false);
							
		inputParams.add(GM_MA);
		inputParams.add(WM_MA);
		inputParams.add(BrainMask_MA);
		inputParams.add(GM_toads);
		inputParams.add(WM_toads);
		inputParams.add(Mask_toads);
		inputParams.add(ErodeMethod);
		inputParams.add(GM_confident_level);

	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		// TODO Auto-generated method stub
		// combined segmentations
		GM_combine=new ParamVolume("combined_GM volume");
		WM_combine=new ParamVolume("combined_WM volume");
		Mask_combine=new ParamVolume("combined_Mask volume");
		
		outputParams.add(GM_combine);
		outputParams.add(WM_combine);
		outputParams.add(Mask_combine);
	}

	@Override
	protected void execute(CalculationMonitor monitor) {
		// TODO Auto-generated method stub

		ImageData GM_MA_vol = GM_MA.getImageData();
		ImageData WM_MA_vol = WM_MA.getImageData();
		ImageData BrainMask_MA_vol = BrainMask_MA.getImageData();
		
		ImageData GM_toads_vol = GM_toads.getImageData();
		ImageData WM_toads_vol = WM_toads.getImageData();
		ImageData Mask_toads_vol = Mask_toads.getImageData();
		
		int rows = GM_MA_vol.getRows();
		int cols = GM_MA_vol.getCols();
		int slices = GM_MA_vol.getSlices();
		
		int i,j,k;
		
		// erode WM_MA
		MASImorph morph = new MASImorph(WM_MA_vol);		
		// get erode image
		if (ErodeMethod.getIndex()==0){
			morph.Minimize3DFilter();
		}else{
			morph.Minimize3DFilterLikeMatlab();
		}
		ImageData WM_MA_vol_erode = morph.GetOutput();
		
		// init output volumes
		ImageDataFloat GM_combine_vol = new ImageDataFloat(rows,cols,slices);
		ImageDataFloat WM_combine_vol = new ImageDataFloat(rows,cols,slices);
		ImageDataInt Mask_combine_vol = new ImageDataInt(rows,cols,slices);
		

		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				for (k = 0; k < slices; k++) {
					if(BrainMask_MA_vol.getInt(i,j,k)>0){ // inside the whole brain mask
						//confident WM
						if(WM_MA_vol_erode.getInt(i,j,k)>0){
							GM_combine_vol.set(i,j,k,0f);
							WM_combine_vol.set(i,j,k,255f);
							Mask_combine_vol.set(i,j,k,255);
						}
						//confident GM
						else if(GM_toads_vol.getFloat(i,j,k)>0
								&&GM_toads_vol.getFloat(i,j,k)<(255*GM_confident_level.getFloat())
								&&GM_MA_vol.getInt(i,j,k)>0){
							GM_combine_vol.set(i,j,k,255f);
							WM_combine_vol.set(i,j,k,0f);
							Mask_combine_vol.set(i,j,k,0);
						}
						//other conditions
						else{
							GM_combine_vol.set(i,j,k,GM_toads_vol.getFloat(i,j,k));
							WM_combine_vol.set(i,j,k,WM_toads_vol.getFloat(i,j,k));
							Mask_combine_vol.set(i,j,k,Mask_toads_vol.getInt(i,j,k));
						}
						
					}else{
						GM_combine_vol.set(i,j,k,0f);
						WM_combine_vol.set(i,j,k,0f);
						Mask_combine_vol.set(i,j,k,0);

					}
					

				}
			}
		}
		
		
		GM_combine_vol.setName("MA_TOADS_GM");
		GM_combine_vol.setHeader(GM_toads_vol.getHeader());
		GM_combine.setValue(GM_combine_vol);
		
		WM_combine_vol.setName("MA_TOADS_WM");
		WM_combine_vol.setHeader(WM_toads_vol.getHeader());
		WM_combine.setValue(WM_combine_vol);
		
		Mask_combine_vol.setName("MA_TOADS_Mask");
		Mask_combine_vol.setHeader(Mask_toads_vol.getHeader());
		Mask_combine.setValue(Mask_combine_vol);


	}

}
