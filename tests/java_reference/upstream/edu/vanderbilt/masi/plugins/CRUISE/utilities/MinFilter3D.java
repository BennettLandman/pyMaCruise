package edu.vanderbilt.masi.plugins.CRUISE.utilities;

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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFileCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

/**
 * @author Yuankai Huo
 * @email yuankai.huo@vanderbilt.edu
 * @version 1.0
 *
 *transfer the toads cruise space image to original space
 *
 */

public class MinFilter3D extends ProcessingAlgorithm{

	private ParamVolume InputVol;
	private ParamVolume OutputVol;
	private ParamOption ErodeMethod;
//	ParamNumberCollection OverlookLabels;
	private ParamString OverlookLabels;
	

	private static final String cvsversion = "$Revision: 1.5 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Conduct imerode like MATLAB" + revnum + "\n";

	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.setPackage("Vanderbilt");
		inputParams.setCategory("format");
		inputParams.setLabel("3D morphological filter");
		inputParams.setName("3D_mParamCollectionorphological_filter");

		// input multi-atlas segmentation
		InputVol=new ParamVolume("Input volume");
		
		// 3D min filter or slice by slice like matlab
		ErodeMethod = new ParamOption("ErodeMethod", new String[] {
				"3D", "3Dmatlab"});
		ErodeMethod.setValue(0);
		ErodeMethod.setMandatory(false);
		
		OverlookLabels = new ParamString("Overlooked Labels in new csf boundary, seperate by comma","120,121,202,203,156,157");
		
		// the labels will not be used in getting csf boundary
//		OverlookLabels = new ParaParamCollectionmNumberCollection("Overlooked Labels in new csf boundary");
//		List<Integer> OverlookList = new ArrayList<Integer>();
//		OverlookList.add(120);
//		OverlookList.add(121);
//		OverlookList.add(202);
//		OverlookList.add(203);
//		OverlookList.add(156);
//		OverlookList.add(157);
//		OverlookLabels.setValue(OverlookList);
//		OverlookLabels.setMandatory(false);
						
		inputParams.add(InputVol);
		inputParams.add(ErodeMethod);
		inputParams.add(OverlookLabels);
		
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		// TODO Auto-generated method stub
		OutputVol = new ParamVolume("Output volume");
		outputParams.add(OutputVol);
	}

	@Override
	protected void execute(CalculationMonitor monitor)
			throws AlgorithmRuntimeException {
		// TODO Auto-generated method stub

		ImageDataInt Inputimg = new ImageDataInt(InputVol.getImageData());
		int rows = Inputimg.getRows();
		int cols = Inputimg.getCols();
		int slices = Inputimg.getSlices();
		
		
		//get labels that will not be used in new csf boundary
//		int []OverlookList;
//		OverlookList = new int[OverlookLabels.size()];
//			for(int i=0;i<OverlookList.length;i++){
//				OverlookList[i] = OverlookLabels.getValue(i).intValue();
//		}

		int[][][] InputVol = Inputimg.toArray3d();
		int[][][] BoundaryVol = new int[rows][cols][slices];
		
		// remove overlooked labels
		InputVol = RemoveLabels(InputVol,OverlookLabels);

		// get unique 133 labels
		Set<Integer> uniqueLabels =  GetUniqueLabels(InputVol);
		Iterator<Integer> itrLabels = uniqueLabels.iterator();		
		while(itrLabels.hasNext()){
			int workingLabel = (int) itrLabels.next();
			if(workingLabel>=100){
				// if the label is >= 100, then get boundary
				BoundaryVol = GetBoundrayForOneLabel(BoundaryVol,InputVol,workingLabel,ErodeMethod.getIndex());
				System.out.println("Done ROI "+workingLabel);// yk add debug
			}
		}
		ImageDataInt Outputimg = new ImageDataInt(BoundaryVol);

		String outputname =  Inputimg.getName()+"_refspace";
		Outputimg.setName(outputname);
		OutputVol.setValue(Outputimg);
		OutputVol.getImageData().setHeader(Inputimg.getHeader());


	}

	public int[][][] RemoveLabels(int[][][]InputVol,ParamString OverlookLabels){
		
		// get removable labels
		String RemoveString = OverlookLabels.getValue();
		String[] RemoveStrings = RemoveString.split(",");
//		int[] RemoveLabels = new int[RemoveStrings.length];
		List<Integer> RemoveLabelsList = new ArrayList<Integer>();
		for (int ii = 0; ii < RemoveStrings.length; ii++) {
		    try {
//		    	RemoveLabels[i] = Integer.parseInt(RemoveStrings[i]);
		    	RemoveLabelsList.add(Integer.parseInt(RemoveStrings[ii]));
		    } catch (NumberFormatException nfe) {};
		}
		
		
		int rows = InputVol.length;
		int cols = InputVol[0].length;
		int slices = InputVol[0][0].length;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					if(RemoveLabelsList.contains(InputVol[i][j][k])){
						InputVol[i][j][k] = 0;
					}
				}
			}
		}
		return InputVol;
	}
	
	
	public int[][][] GetBoundrayForOneLabel(int[][][] BoundaryVol,int[][][] InputVol,int workingLabel,int ErodeMethod){
		int rows = InputVol.length;
		int cols = InputVol[0].length;
		int slices = InputVol[0][0].length;
		int[][][] ErodeInputVol;

		int[][][] BinaryInputVol = new int[rows][cols][slices];
		// intialize the  BinaryInputVol
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					if(InputVol[i][j][k] == workingLabel){
						BinaryInputVol[i][j][k] = 1;
					}else{
						BinaryInputVol[i][j][k] = 0;
					}

				}
			}
		}

		// get erode image
		if (ErodeMethod==0){
			ErodeInputVol = Minimize3DFilter(BinaryInputVol);
		}else{
			ErodeInputVol = Minimize3DFilterLikeMatlab(BinaryInputVol);
		}

		// get boundary image by substract erode image from original image
		for (int i = 1; i < rows-1; i++) {
			for (int j = 1; j < cols-1; j++) {
				for (int k = 1; k < slices-1; k++) {
					if (BinaryInputVol[i][j][k]!=ErodeInputVol[i][j][k]){
						BoundaryVol[i][j][k] = getTrueBoundaryValue(InputVol,i,j,k);
//						BoundaryVol[i][j][k] = 1;
					}

				}
			}
		}


		return BoundaryVol;
	}

	public int getTrueBoundaryValue(int[][][] InputVol,int i,int j,int k){
		int TrueBoundaryValue;
		int x,y,z;
		int centerVoxel = InputVol[i][j][k];
		boolean finddifferentVoxel = false;  //true means find other labels other than center Voxels

		for (x = i-1; x <= i+1; x++) {
			for (y = j-1; y <= j+1; y++) {
				for (z = k-1; z <= k+1; z++) {
					if (InputVol[x][y][z] != 0){
						if (InputVol[x][y][z] < 100){
							TrueBoundaryValue = 0;
							return TrueBoundaryValue; // touch subcortical regions, not good
						}

						if (InputVol[x][y][z]!=centerVoxel){
							finddifferentVoxel = true;
						}
					}
				}
			}
		}

		if(finddifferentVoxel){
			TrueBoundaryValue = 1;
		}else{
			TrueBoundaryValue = 0;
		}
		return TrueBoundaryValue;
	}

	public Set<Integer> GetUniqueLabels(int[][][] inputVol){
		int rows = inputVol.length;
		int cols = inputVol[0].length;
		int slices = inputVol[0][0].length;
		Set<Integer> uniques = new HashSet<Integer>();

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					uniques.add(inputVol[i][j][k]);
				}
			}
		}

		return uniques;
	}



	public ImageData Minimize3DFilter(ImageData Inputimg){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int rows = Inputimg.getRows();
		int cols = Inputimg.getCols();
		int slices = Inputimg.getSlices();
		int i,j,k;
		float PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		float TempMin;

		ImageData outputimg = Inputimg.clone();

		//		float tempval;

		for (i = 1; i < rows-1; i++) {
			for (j = 1; j < cols-1; j++) {
				for (k = 1; k < slices-1; k++) {

					// disk 3D mask
					PointOrigin = Inputimg.getFloat(i, j, k);
					PointNorth = Inputimg.getFloat(i-1, j, k);
					PointSouth = Inputimg.getFloat(i+1, j, k);
					PointWest = Inputimg.getFloat(i, j-1, k);
					PointEast = Inputimg.getFloat(i, j+1, k);
					PointFront = Inputimg.getFloat(i, j, k-1);
					PointBack = Inputimg.getFloat(i, j, k+1);

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.min(TempMin,PointNorth);	
					TempMin = Math.min(TempMin,PointSouth);	
					TempMin = Math.min(TempMin,PointWest);	
					TempMin = Math.min(TempMin,PointEast);	
					TempMin = Math.min(TempMin,PointFront);	
					TempMin = Math.min(TempMin,PointBack);	

					outputimg.set(i, j, k, TempMin);
				}
			}
		}

		return outputimg;
	}

	public int[][][] Minimize3DFilterLikeMatlab(int[][][] InputVol){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int rows = InputVol.length;
		int cols = InputVol[0].length;
		int slices = InputVol[0][0].length;
		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		int[][][] outputimg = new int[rows][cols][slices];

		//		float tempval;
		for (k = 1; k < slices-1; k++) {
			for (i = 1; i < rows-1; i++) {
				for (j = 1; j < cols-1; j++) {

					// disk 3D mask
					PointOrigin = InputVol[i][j][k];
					PointNorth = InputVol[i-1][j][k];
					PointSouth = InputVol[i+1][j][k];
					PointWest = InputVol[i][j-1][k];
					PointEast = InputVol[i][j+1][k];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.min(TempMin,PointNorth);	
					TempMin = Math.min(TempMin,PointSouth);	
					TempMin = Math.min(TempMin,PointWest);	
					TempMin = Math.min(TempMin,PointEast);	


					outputimg[i][j][k] = TempMin;
				}
			}
		}

		return outputimg;
	}

	public int[][][] Minimize3DFilter(int[][][] InputVol){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int rows = InputVol.length;
		int cols = InputVol[0].length;
		int slices = InputVol[0][0].length;
		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		int[][][] outputimg = new int[rows][cols][slices];

		//		float tempval;

		for (i = 1; i < rows-1; i++) {
			for (j = 1; j < cols-1; j++) {
				for (k = 1; k < slices-1; k++) {

					// disk 3D mask
					PointOrigin = InputVol[i][j][k];
					PointNorth = InputVol[i-1][j][k];
					PointSouth = InputVol[i+1][j][k];
					PointWest = InputVol[i][j-1][k];
					PointEast = InputVol[i][j+1][k];
					PointFront = InputVol[i][j][k-1];
					PointBack = InputVol[i][j][k+1];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.min(TempMin,PointNorth);	
					TempMin = Math.min(TempMin,PointSouth);	
					TempMin = Math.min(TempMin,PointWest);	
					TempMin = Math.min(TempMin,PointEast);	
					TempMin = Math.min(TempMin,PointFront);	
					TempMin = Math.min(TempMin,PointBack);	

					outputimg[i][j][k] = TempMin;
				}
			}
		}

		return outputimg;
	}

}
