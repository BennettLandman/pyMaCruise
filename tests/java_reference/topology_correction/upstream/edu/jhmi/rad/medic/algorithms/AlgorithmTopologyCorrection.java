package edu.jhmi.rad.medic.algorithms;
import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.structures.BinaryTree;
import edu.jhmi.rad.medic.utilities.*;

import gov.nih.mipav.model.algorithms.AlgorithmBase;
import gov.nih.mipav.model.structures.*;
import gov.nih.mipav.view.*;

import java.io.*;
import java.util.*;
import java.lang.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import WildMagic.LibFoundation.Mathematics.Vector3f;


/**
 *
 *  algorithm for correcting the topology of binary and scalar images.
 *  <p>
 *	It is based on a fast marching method with topology-preserving constraints.
 *	This version is publicly released and stable.
 *
 *	@version    April 2005
 *	@author     Pierre-Louis Bazin
 *  @see 		JDialogTopologyCorrection
 *  @see 		TopologyPropagation
 *		
 *
*/
public class AlgorithmTopologyCorrection extends AlgorithmBase {
	private static final String cvsversion = "$Revision: 1.1 $";
	public static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");

	public static String get_version() {
		return revnum;
	}


    // Fuzzy images require 1 image for each class
    // Hard images 1 image with assigned clusters
    private ModelImage  destImage;
    private ViewUserInterface       userInterface;

    
     // information on the parts of the image to process
	private boolean 	useWholeImage = true;
	private boolean 	cropBackground = false;
    private float       signalThreshold=0.0f;
    private int         nxmin,nymin,nzmin,nxmax,nymax,nzmax;

    // image size
	private int 		nx,ny,nz; // original image dimensions
	private int			mx,my,mz; // reduced image dimensions
    private int         x0,xN,y0,yN,z0,zN;
	private float		Imin, Imax; // original image extents
	
    // algorithm parameters
    private float   	minDistance;
	private int         dim;
    private int         cObj;
    private int         cBack;
    private String      type;
    private	String		inputSkel;
    private BitSet      inputPaint;
	private float   	highestLevel = 1.0f;
    private float   	lowestLevel  = 0.0f;
	
	// other constants
    private float   	maxDistance = 0.1f;
    private boolean     useMinMax = true;
    private boolean			extra = false;
    private ModelImage		extraImage;
	private	String		propagDir ="downward";
	private Vector3f[]  inputSeed = null;
    private	boolean		thresholdStop = false;
	private int         maxObj = 100;
    
	private boolean		verbose = false;
	
    /**
    *	Constructor for 3D images in which changes are placed in a predetermined destination image.
    *   @param destImage_      Image model where result image is to stored.
    *   @param srcImage_       Source image model.
    */
	public AlgorithmTopologyCorrection(ModelImage destImage_, ModelImage srcImage_,
                                    String type_,
                                    float minDist_, 
									int dim_, int cO_, int cB_, 
									String skel_, BitSet paint_,
									float low_, float high_, boolean minmax_,
									String propagType_) {
        super(null, srcImage_);
		userInterface = ViewUserInterface.getReference();
        destImage = destImage_;  // Put results in destination image.
        
        type = type_;
		minDistance = minDist_;
        dim = dim_;
        cObj = cO_;
        cBack = cB_;
        inputSkel = skel_;
        inputPaint = paint_;
		lowestLevel = low_;
		highestLevel = high_;
		useMinMax = minmax_;
		
		if ( (propagType_.equals("background->object")) && (inputSkel.equals("intensity")) )
			propagDir = "upward";
		else
			propagDir = "downward";
		
	}


    /**
    *	Prepares this class for destruction.
    */
	public void finalize(){
	    destImage   = null;
	    srcImage    = null;
        System.gc();
        super.finalize();
	}

    /**
    *	Constructs a string of the contruction parameters and outputs the string to the messsage frame if the logging
    *   procedure is turned on.
    *	@param destinationFlag	If true the log includes the name of the destination flag.
    */
	/*
    private void constructLog(boolean destinationFlag) {
        if ( destinationFlag == false) {
            historyString = new String( "TopologyCorrection(" + ")");
        }
        else  {
            historyString = new String( "TopologyCorrection(" + ")");
        }
        historyString += "\n";  // append a newline onto the string
        writeLog();
    }
	*/

    /**
     * Writes the logString to the appropriate log area. Overrides the AlgorithmBase <code>writeLog()</code> to append
     * the history to all of the destination images we've created.
     */
	 /*
    protected void writeLog() {
        // write to the history area
        if (Preferences.is(Preferences.PREF_LOGGING_ENABLED) && isCompleted()) {
			if ((destImage != null) && (destImage.getHistoryArea() != null)) {
				if (srcImage != null) {
					destImage.getHistoryArea().setText(srcImage.getHistoryArea().getText());
                }
				if (historyString != null) {
                    destImage.getHistoryArea().append(historyString);
				}
            } else if (srcImage != null) {
				if (historyString != null) {
                    srcImage.getHistoryArea().append(historyString);
                }
            }
        }
    }
	*/
	
	/**
    *   Starts the algorithm.
    */
	public void runAlgorithm() {

        if (srcImage  == null) {
            displayError("Source Image is null");
            //notifyListeners(this);
            return;
        }
        if (destImage  == null) {
            displayError("Source Image is null");
            //notifyListeners(this);
            return;
        }

        // start the timer to compute the elapsed time
        setStartTime();

        if (destImage != null){     // if there exists a destination image
            //constructLog(true);
			
			float [] buffer;
            if (srcImage.getNDims() == 2){
				try {
					// image length is length in 2 dims
					int length = srcImage.getExtents()[0]  * srcImage.getExtents()[1];
					buffer       = new float[length];
					srcImage.exportData(0,length, buffer); // locks and releases lock
					fireProgressStateChanged(srcImage.getImageName(), "Processing image ...");
				}
				catch (IOException error) {
					buffer = null;
					errorCleanUp("Algorithm Topology Correction reports: source image locked", true);
					return;
				}
				catch (OutOfMemoryError e){
					buffer = null;
					errorCleanUp("Algorithm Topology Correction reports: out of memory", true);
					return;
				}

				// init dimensions
				nx = srcImage.getExtents()[0];
				ny = srcImage.getExtents()[1];
				nz = 1;
				
				// main algorithm
				calcTopology(buffer);
				
			}
	        else if (srcImage.getNDims() > 2) {
			   try {
					// image length is length in 3 dims
					int length = srcImage.getExtents()[0] *
							 srcImage.getExtents()[1] *
							 srcImage.getExtents()[2];
					buffer = new float[length];
					srcImage.exportData(0,length, buffer); // locks and releases lock
					fireProgressStateChanged(srcImage.getImageName(), "Processing image ...");
				}
				catch (IOException error) {
					buffer = null;
					errorCleanUp("Algorithm Topology Correction: source image locked", true);
					return;
				}
				catch (OutOfMemoryError e){
					buffer = null;
					errorCleanUp("Algorithm Topology Correction: Out of memory creating process buffer", true);
					return;
				}
		
				// init dimensions
				nx = srcImage.getExtents()[0];
				ny = srcImage.getExtents()[1];
				nz = srcImage.getExtents()[2];
				
				// main algorithm
				calcTopology(buffer);
				
			}
        }

        // compute the elapsed time
        computeElapsedTime();
    } // end run()

    /**
    *	computes the topology correction for 3D images
    */
    private void calcTopology(float img[]){
		float[][][]   	result=null; 
        float[][][]     image;
		boolean[][][]     objectMask;
		boolean[][][]     obj = null;
		boolean[][]     obj2D = null;
		int[][][]     label = null;
		int[][]     label2D = null;
		int 		x,y,z,k,m,i,j,l;
		float   	num,den;
		float   	ngb, neighbors;
		float   	dist;
        float[][][]   	ext = null; 
        float[]     buffer;
		boolean[][][]     paint;
		TopologyPropagation algorithm = null;
		String info;
		float[][][] lb = null;
		float level,Dmax;
		int maxLabels,Nc;
		float[] labels;
		boolean isNew;
		float best;
        int xM,yM,zM;
		boolean stop;
                
			            
		int mod;
		int progress;

		long start_time, inner_loop_time;

		if (verbose) MedicUtilPublic.displayMessage("start correction\n");
				
		if (verbose) MedicUtilPublic.displayMessage("extract data..\n");
		
		/* pre-processing : expand boundaries, so that we don't have to worry for them */
        try {
            expandSize();
            image = expandBoundaries(img);
            if (inputSkel.equals("paint_mask")) paint = expandBoundaries(inputPaint);
            else paint = null;
            // dispose of img
            img = null;
        } catch (OutOfMemoryError e) {
            img = null;
            errorCleanUp("Algorithm Topology Correction: Out of memory creating process buffer", true);
            setCompleted(false);
            return;
        }
		if (verbose) MedicUtilPublic.displayMessage("set parameters..\n");

		// compute exact min and max (calcMinMax() returns integers)
		Imin = 1e30f;
		Imax = -1e30f;
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			if (image[x][y][z] < Imin) Imin = image[x][y][z];
			if (image[x][y][z] > Imax) Imax = image[x][y][z];
		}
		if (verbose) MedicUtilPublic.displayMessage("min: "+Imin+", max:"+Imax+"\n");
		
		// use relative values for the levels
		if (useMinMax) {
			highestLevel = Imax - (1.0f-highestLevel)*(Imax-Imin);
			lowestLevel  = Imin + lowestLevel*(Imax-Imin);
		}
		if (verbose) MedicUtilPublic.displayMessage("high: "+highestLevel+", low:"+lowestLevel+"\n");
		
		if (verbose) MedicUtilPublic.displayMessage("create mask..\n");
		
		// process only parts where intensity > lowest 
		// (shrink the image, but doesn't use the mask for topology purposes)
		objectMask = createObjectMask(image, lowestLevel);
        
		// create smaller images for speed and memory
		if (verbose) MedicUtilPublic.displayMessage("compute boundaries..\n");
		computeProcessingBoundaries(objectMask);
		
		if (verbose) MedicUtilPublic.displayMessage("shrink images..\n");
		image = reduceImageSize(image); 
		objectMask = reduceImageSize(objectMask);
        if (inputSkel.equals("paint_mask")) paint = reduceImageSize(paint);
        
		if (verbose) MedicUtilPublic.displayMessage("new dimensions: "+mx+"x"+my+"x"+mz+"\n");
		
		/////////////////////////////////////////
		//          MAIN ALGORITHM             //
		/////////////////////////////////////////
		if (verbose) MedicUtilPublic.displayMessage("start algorithm\n");
 
		// record time
		 start_time = System.currentTimeMillis();
		        
		fireProgressStateChanged("topology correction");
		System.out.println("Start ALGORITHM");
		// correction of the topology for a scalar image			
        if (type.equals("scalar_image")) {
			// use relative values for the levels
			if (useMinMax) {
				minDistance = minDistance*(Imax-Imin);
			}
			// look for the starting point
			if (inputSkel.equals("intensity")) {
				if (propagDir.equals("downward")) {
					// create a label image for the object (value >= middleLevel)
					obj = ObjectProcessing.objectFromImage(image, mx, my, mz,
															0.5f*(lowestLevel+highestLevel), 
															highestLevel, ObjectProcessing.SUPEQUAL, ObjectProcessing.INFERIOR);
					
					if ((cObj==6) || (cObj==4)) label = ObjectProcessing.connected6Object3D(obj, mx, my, mz);
					else if ((cObj==18) || (cObj==8)) label = ObjectProcessing.connected18Object3D(obj, mx, my, mz);
					else if (cObj==26)  label = ObjectProcessing.connected26Object3D(obj, mx, my, mz);
					
					// find the largest label
					obj = null;
					obj = ObjectProcessing.largestObjectFromLabel(label, mx,my,mz);
					// set the starting point for propagation
					best = lowestLevel;
					xM=0;yM=0;zM=0;
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						if ( (obj[x][y][z]) && (image[x][y][z]>best) ) {
							best = image[x][y][z];
							xM=x; yM=y; zM=z;
						}
					}
					// create a paint mask with just this point
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					paint[xM][yM][zM] = true;
					if (verbose) MedicUtilPublic.displayMessage("starting point: "+image[xM][yM][zM]+"("+xM+","+yM+","+zM+"\n");
				} else if (propagDir.equals("upward")) {
					// create a paint mask with the outer boundary
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					for (x=0;x<mx;x++) for (y=0;y<my;y++) {
						paint[x][y][0] = true;
						paint[x][y][1] = true;
						paint[x][y][mz-1] = true;
						paint[x][y][mz-2] = true;
					}
					for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[0][y][z] = true;
						paint[1][y][z] = true;
						paint[mx-1][y][z] = true;
						paint[mx-2][y][z] = true;
					}
					for (x=0;x<mx;x++) for (z=0;z<mz;z++) {
						paint[x][0][z] = true;
						paint[x][1][z] = true;
						paint[x][my-1][z] = true;
						paint[x][my-2][z] = true;
					}
					if (verbose) MedicUtilPublic.displayMessage("starting from a bounding box\n");
				}
            } else if (inputSkel.equals("paint_mask")) {
				// set all points in the mask to highestLevel
				for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
					if (paint[x][y][z]) image[x][y][z] = highestLevel;
				}
			}
			// flatten the values above and below the levels
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z] > highestLevel) image[x][y][z] = lowestLevel;
				if (image[x][y][z] < lowestLevel)  image[x][y][z] = lowestLevel;
			}
			System.out.println("TOPOLOGY PROPOGATION");
			// start the algorithm
			algorithm = new TopologyPropagation(image, objectMask, mx, my, mz, lowestLevel, highestLevel, 
												0.0f, minDistance, dim, cObj, cBack, true, 
												"paint mask", null, paint, userInterface, getProgressChangeListener());
			System.out.println("FINISHED PROPOGATION");
			if (verbose) MedicUtilPublic.displayMessage("initialisation ("+algorithm.isWorking()+")\n");
												
			// initialize from the paint mask
			// redundant ? algorithm.initDownwardPaintLabels();
			
			if (propagDir.equals("upward")) 
				algorithm.propagateUpwardExactSmoothing();
			else 
				algorithm.propagateDownwardExactSmoothing();
			
			result = extendImageSize(algorithm.exportDistance());
		
			algorithm.finalize();
			
		} else if (type.equals("binary_object")) {
			// extract the distance function or the object
			obj = ObjectProcessing.objectFromImage(image,mx,my,mz,lowestLevel,highestLevel,ObjectProcessing.SUPERIOR,ObjectProcessing.INFEQUAL);
			ext = new float[mx][my][mz];
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (obj[x][y][z]) ext[x][y][z] = 1.0f;
				else ext[x][y][z] = 0.0f;
			}
			// start the algorithm	
			if (propagDir.equals("upward")) {
				algorithm = new TopologyPropagation(ext, objectMask, mx, my, mz, 
											0.0f, 0.5f, 0.0f, 0.0f,
											dim, cObj, cBack, false, 
											"intensity", null, null, userInterface, getProgressChangeListener());
	 
				// compute distance function
				algorithm.propagateDownwardUnconstrainedGeometricDistance();
				image = algorithm.exportDistance();
				algorithm.finalize();			
			} else {
				algorithm = new TopologyPropagation(ext, objectMask, mx, my, mz, 
											0.5f, 1.0f, 0.0f, 0.0f,
											dim, cObj, cBack, false, 
											"intensity", null, null, userInterface, getProgressChangeListener());
	 
				// compute distance function
				algorithm.propagateUpwardUnconstrainedGeometricDistance();
				image = algorithm.exportDistance();
				algorithm.finalize();
			}
			// reset highest and lowest to fit the distance function
			highestLevel = 0.0f;
			lowestLevel = 0.0f;
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z]>highestLevel) highestLevel=image[x][y][z];
				if (image[x][y][z]<lowestLevel)  lowestLevel=image[x][y][z];
			}		
			// use relative values for the levels
			if (useMinMax) {
				minDistance = minDistance*(Imax-Imin);
			}
					
			// find the highest point to start from
            if (inputSkel.equals("intensity")) {
				if (propagDir.equals("downward")) {
					// set the starting point for propagation
					best = 0.0f;
					xM=0;yM=0;zM=0;
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						if (image[x][y][z]>best) {
							best = image[x][y][z];
							xM=x; yM=y; zM=z;
						}
					}
					// create a paint mask with just this point
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					paint[xM][yM][zM] = true;
					if (verbose) MedicUtilPublic.displayMessage("starting point: "+image[xM][yM][zM]+"("+xM+","+yM+","+zM+"\n");
				} else if (propagDir.equals("upward")) {
					// create a paint mask with the outer boundary
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					for (x=0;x<mx;x++) for (y=0;y<my;y++) {
						paint[x][y][0] = true;
						paint[x][y][1] = true;
						paint[x][y][mz-1] = true;
						paint[x][y][mz-2] = true;
					}
					for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[0][y][z] = true;
						paint[1][y][z] = true;
						paint[mx-1][y][z] = true;
						paint[mx-2][y][z] = true;
					}
					for (z=0;z<mz;z++) for (x=0;x<mx;x++) {
						paint[x][0][z] = true;
						paint[x][1][z] = true;
						paint[x][my-1][z] = true;
						paint[x][my-2][z] = true;
					}
					if (verbose) MedicUtilPublic.displayMessage("starting from a bounding box\n");
										
				}
            } else if (inputSkel.equals("paint_mask")) {
				// set all points in the mask to highestLevel
				for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
					if (paint[x][y][z]) image[x][y][z] = highestLevel;
				}
			}
			// flatten the values above and below the levels
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z] > highestLevel) image[x][y][z] = lowestLevel;
				if (image[x][y][z] < lowestLevel)  image[x][y][z] = lowestLevel;
			}
			
			// start the algorithm
			algorithm = new TopologyPropagation(image, objectMask, mx, my, mz, lowestLevel, highestLevel, 
												0.0f, minDistance, dim, cObj, cBack, true, 
												"paint mask", null, paint, userInterface, getProgressChangeListener());

			if (verbose) MedicUtilPublic.displayMessage("initialisation ("+algorithm.isWorking()+")\n");
												
			// initialize from the paint mask
			// redundant ? algorithm.initDownwardPaintLabels();
			
			if (propagDir.equals("upward")) 
				algorithm.propagateUpwardExactSmoothing();
			else 
				algorithm.propagateDownwardExactSmoothing();
			
			result = extendImageSize(algorithm.exportDistance());
		
			algorithm.finalize();	
		}					

        if (threadStopped) {
			algorithm.finalize();
			finalize();
			setCompleted(false);
        	return;
        } 

		// debug
		if (verbose) MedicUtilPublic.displayMessage("total time: (milliseconds): " + (System.currentTimeMillis()-start_time)); 

		// extract results (segmentation and/or classes) and store them in destImage[]
		fireProgressStateChanged("creating result images...");
		nx = nx-2;
		ny = ny-2;
		nz = nz-2;
		try {
			buffer = new float[nx*ny*nz];
			for (x=0;x<nx;x++) {
				for (y=0;y<ny;y++) {
					for (z=0;z<nz;z++) {
						buffer[x + nx*y + nx*ny*z] = result[x+1][y+1][z+1];
					}
				}
			}
			destImage.importData(0, buffer, true);
		} catch (OutOfMemoryError e) {
			buffer = null;
			errorCleanUp("Algorithm Topology: Out of memory creating hard classification", true);
			finalize();
			setCompleted(false);
			return;
		} catch (IOException error) {
			errorCleanUp("Algorithm Topology: export problem to destImage[]", true);
			finalize();
			setCompleted(false);
			return;
		}

        image = null;
		objectMask = null;
		result = null;
		lb = null;
		obj = null;
		label = null;
		System.out.println("COMPLETED");
		setCompleted(true);
    } // calcTopology3D
	
	/** expand boundaries for spatial comutations */
	private void expandSize() {
		nx = nx+2;
		ny = ny+2;
		nz = nz+2;
	}			
	private float[][][] expandBoundaries(float[] image) {
		int 		x,y,z;
		float[][][] 	tmp;
		
		tmp = new float[nx][ny][nz];
		for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++)
					tmp[x][y][z] = image[ (x-1) + (nx-2)*(y-1) + (nx-2)*(ny-2)*(z-1) ];
		
		return tmp;
	}
	private boolean[][][] expandBoundaries(BitSet image) {
		int 		x,y,z;
		boolean[][][] 	tmp;
		
		tmp = new boolean[nx][ny][nz];
		for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++)
					tmp[x][y][z] = image.get( (x-1) + (nx-2)*(y-1) + (nx-2)*(ny-2)*(z-1) );
		
		return tmp;
	}
	
	/** creates a mask for unused data */
	private boolean[][][] createObjectMask(float[][][] image, float val) {
		int 		x,y,z;
		boolean[][][]  	objMask;

		// uses only values over the threshold, if mask used
		objMask = new boolean[nx][ny][nz];
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++)
				for (z=0;z<nz;z++) {
                    if (useWholeImage) {
                        if ( (cropBackground) && (image[x][y][z] <= val) )
                            objMask[x][y][z] = false;
                        else
                            objMask[x][y][z] = true;
                    } else if (mask.get(x+nx*y+nx*ny*z) ) {
                        if ( (cropBackground) && (image[x][y][z] <= val) )
                            objMask[x][y][z] = false;
                        else
                            objMask[x][y][z] = true;
                    }
                }
		// remove the boundary from the computations
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++) {
				objMask[x][y][0] = false;
				objMask[x][y][nz-1] = false;
			}
		for (y=0;y<ny;y++)
			for (z=0;z<nz;z++) {
				objMask[0][y][z] = false;
				objMask[nx-1][y][z] = false;
			}
		for (z=0;z<nz;z++)
			for (x=0;x<nx;x++) {
				objMask[x][0][z] = false;
				objMask[x][ny-1][z] = false;
			}

		return objMask;
	} // createObjectMask
        
    /** sets the processing lower and upper boundaries */
    private void computeProcessingBoundaries(boolean[][][] objMask) {
		int 		x,y,z;
        
        x0 = nx;
        xN = 0;
        y0 = ny;
        yN = 0;
        z0 = nz;
        zN = 0;
        for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++) {
                    if (objMask[x][y][z]) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
		// update the smaller size parameters (include am extended boundary)
		mx = xN - x0 + 1 + 2;
		my = yN - y0 + 1 + 2;
		mz = zN - z0 + 1 + 2;
		
        // debug
        System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
	/** create smaller image (for saving memory) */
	private float[][][] reduceImageSize(float[][][] image) {
		float[][][] smaller = new float[mx][my][mz];

		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					smaller[x-x0+1][y-y0+1][z-z0+1] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
	private boolean[][][] reduceImageSize(boolean[][][] image) {
		boolean[][][] smaller = new boolean[mx][my][mz];

		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					smaller[x-x0+1][y-y0+1][z-z0+1] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
	private void reducePointPosition(Vector3f[] pt) {
		for (int n=0;n<pt.length;n++) {
            // from image to image+boundary to reduced image
            pt[n].X = pt[n].X - x0 + 1 + 2;
            pt[n].Y = pt[n].Y - y0 + 1 + 2;
            pt[n].Z = pt[n].Z - z0 + 1 + 2;
		}
		return;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(float[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(int[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = (float)image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(byte[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = (float)image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
    /** retrieve original size from smaller image */
	private float[][][] extendImageSize(boolean[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
                    if (image[x-x0+1][y-y0+1][z-z0+1])
                        larger[x][y][z] = 1.0f;
                    else
                        larger[x][y][z] = 0.0f;
				}
			}
		}
		return larger;
	}
}
