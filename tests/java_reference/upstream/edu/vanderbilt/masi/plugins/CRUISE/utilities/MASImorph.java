package edu.vanderbilt.masi.plugins.CRUISE.utilities;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;

/**
 * @author Yuankai Huo
 * @email yuankai.huo@vanderbilt.edu
 * @version 1.0
 *
 *morphological operation
 *1. emrode
 *
 */


public class MASImorph {

	protected int vol3d[][][] = null;
	protected ImageHeader vol3d_header = null;

	protected ImageData output_vol = null;
	protected int output_vol3d [][][] =null;
	private  int rows,cols,slices;

	public MASImorph(ImageData vol){
		ImageDataInt Inputimg = new ImageDataInt(vol);
		vol3d = Inputimg.toArray3d();
		vol3d_header = Inputimg.getHeader();
		// rows, cols, slices
		rows = Inputimg.getRows();
		cols = Inputimg.getCols();
		slices = Inputimg.getSlices();;
		// init output
		output_vol3d = new int[rows][cols][slices];
	}

	public MASImorph(int[][][] voldata){
		vol3d = voldata;
		// rows, cols, slices
		rows = voldata.length;
		cols = voldata[0].length;
		slices = voldata[0][0].length;
	}

	public ImageData GetOutput(){
		ImageData output_vol = new ImageDataInt(output_vol3d); 
		output_vol.setHeader(vol3d_header);
		return output_vol;
	}

	public void Minimize3DFilterLikeMatlab(){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		//		float tempval;
		for (k = 1; k < slices-1; k++) {
			for (i = 1; i < rows-1; i++) {
				for (j = 1; j < cols-1; j++) {

					// disk 3D mask
					PointOrigin = vol3d[i][j][k];
					PointNorth = vol3d[i-1][j][k];
					PointSouth = vol3d[i+1][j][k];
					PointWest = vol3d[i][j-1][k];
					PointEast = vol3d[i][j+1][k];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.min(TempMin,PointNorth);	
					TempMin = Math.min(TempMin,PointSouth);	
					TempMin = Math.min(TempMin,PointWest);	
					TempMin = Math.min(TempMin,PointEast);	


					output_vol3d[i][j][k] = TempMin;
				}
			}
		}

	}

	public void Minimize3DFilter(){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1

		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		//		float tempval;

		for (i = 1; i < rows-1; i++) {
			for (j = 1; j < cols-1; j++) {
				for (k = 1; k < slices-1; k++) {

					// disk 3D mask
					PointOrigin = vol3d[i][j][k];
					PointNorth = vol3d[i-1][j][k];
					PointSouth = vol3d[i+1][j][k];
					PointWest = vol3d[i][j-1][k];
					PointEast = vol3d[i][j+1][k];
					PointFront = vol3d[i][j][k-1];
					PointBack = vol3d[i][j][k+1];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.min(TempMin,PointNorth);	
					TempMin = Math.min(TempMin,PointSouth);	
					TempMin = Math.min(TempMin,PointWest);	
					TempMin = Math.min(TempMin,PointEast);	
					TempMin = Math.min(TempMin,PointFront);	
					TempMin = Math.min(TempMin,PointBack);	

					output_vol3d[i][j][k] = TempMin;
				}
			}
		}

	}

	public void Maximum3DFilterLikeMatlab(){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1
		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		//		float tempval;
		for (k = 1; k < slices-1; k++) {
			for (i = 1; i < rows-1; i++) {
				for (j = 1; j < cols-1; j++) {

					// disk 3D mask
					PointOrigin = vol3d[i][j][k];
					PointNorth = vol3d[i-1][j][k];
					PointSouth = vol3d[i+1][j][k];
					PointWest = vol3d[i][j-1][k];
					PointEast = vol3d[i][j+1][k];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.max(TempMin,PointNorth);	
					TempMin = Math.max(TempMin,PointSouth);	
					TempMin = Math.max(TempMin,PointWest);	
					TempMin = Math.max(TempMin,PointEast);	


					output_vol3d[i][j][k] = TempMin;
				}
			}
		}

	}
	
	public void Maximum3DFilter(){
		//fMethod : 1 means add with cap one, 2 means with cap 0.1

		int i,j,k;
		int PointOrigin,PointNorth,PointSouth,PointWest,PointEast,PointFront,PointBack;
		int TempMin;

		//		float tempval;

		for (i = 1; i < rows-1; i++) {
			for (j = 1; j < cols-1; j++) {
				for (k = 1; k < slices-1; k++) {

					// disk 3D mask
					PointOrigin = vol3d[i][j][k];
					PointNorth = vol3d[i-1][j][k];
					PointSouth = vol3d[i+1][j][k];
					PointWest = vol3d[i][j-1][k];
					PointEast = vol3d[i][j+1][k];
					PointFront = vol3d[i][j][k-1];
					PointBack = vol3d[i][j][k+1];

					//get min value
					TempMin = PointOrigin;	
					TempMin = Math.max(TempMin,PointNorth);	
					TempMin = Math.max(TempMin,PointSouth);	
					TempMin = Math.max(TempMin,PointWest);	
					TempMin = Math.max(TempMin,PointEast);	
					TempMin = Math.max(TempMin,PointFront);	
					TempMin = Math.max(TempMin,PointBack);	

					output_vol3d[i][j][k] = TempMin;
				}
			}
		}

	}
}
