package edu.jhmi.rad.medic.libraries;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.*;
import edu.jhmi.rad.medic.utilities.*;
import WildMagic.LibFoundation.Mathematics.*;

import edu.jhmi.rad.medic.structures.BinaryTree;
 
/**
 *
 *  This algorithm propagates topology in images.
 *  <p>
 *  Different options are possible:
 *  approximate or exact, upward or downward smoothing.
 *  The propagation is based on a fast marching technique, and applies 
 *  to distance functions or grayscale images.
 *	
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class TopologyPropagation {
	
	// propagation flags, also used for checking topology
	// ! (SKELETON,OBJECT,TEST) must be > others
	public  static final byte    BACKFRONTIER = 110;
	private static final byte    BACKMASK = 	100;
	public  static final byte    OBJECT =       80;
	private static final byte    TEST =         70;
	private static final byte    BOUNDARY =     40;
	private static final byte    BACKGROUND =   30;
	private static final byte    CRITICAL =     20;
	private static final byte    MBACK =        18;
	private static final byte    MBOUND =       16;
	private static final byte    MTEST =        14;
	private static final byte    MOBJ =         12;
	private static final byte    MASK =         10;
	private static final byte    FRONTIER =     00;
	
	private static final byte	 OBJSET		=	8;
	private static final byte	 OBJTEST	=	7;
	private static final byte	 OBJBOUND	=	6;
	private static final byte	 OBJCRIT	=	5;
	private static final byte	 OBJMASK	=	3;
	
	private static final byte	 UNSET		=	4;
	
	private static final byte	 BGMASK		=	5;
	private static final byte	 BGCRIT		=	3;
	private static final byte	 BGBOUND	=	2;
	private static final byte	 BGTEST		=	1;
	private static final byte	 BGSET		=	0;
	
	
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-6f;
	
	// data buffers
	private 	float[][][] 	image;  			// original image
	private 	float[][][]     dist;   			// propagated function
	private 	byte[][][]      label;  			// labels for the marching progression
	private 	boolean[][][]   mask;   			// image mask: true for data points
	private 	float[][][]     geom;   			// a geometric distance term for regularisation
	private 	byte[][][]		critical;   		// a map of the critical points
	private 	BinaryTree  	boundary;   		// the binary used for the fast marching
	private static	int 		nx,ny,nz;   		// image dimensions
	
	// parameters
	private 	float   		upLevel;            // the levelset to start the upward algorithms
	private 	float   		downLevel;          // the levelset to start the downward algorithms
	private		int				dimension = 2;      // give the dimension of the computation, 2D or 3D
    private     int             connectObj = 6;     // the connectivity of the foreground object
    private     int             connectBack = 26;   // the connectivity of the background
	private		float			maxDistance;		// the maximum change allowed from an original image
	private		boolean			thresholdStop = false;	// if we use the threshold to stop the propagation
	private     String          inputSkel;          // the type of input for the skeleton: intensity, seed point or seed region
    private     Vector3f[]      inputSeed;          // the seed points for seed skeleton
    private     boolean[][][]   inputPaint;         // the mask for paint skeleton
	private		float			minDistance;		// the minimum distance allowed between neighboring points
	private		int			diagnosticQueueEvents = 0;
    
	// computation variables
	private		int				obj;
	private		int				back;
	private		int				dz;
	
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	ViewUserInterface			UI;
    ViewJProgressBar            progressBar;
	 boolean		debug=false;
	 boolean		verbose=true;
	/**
	 *  constructor
	 */
	public TopologyPropagation(float[][][] image_, boolean [][][] mask_, 
									int nx_, int ny_, int nz_,
									float upLevel_, float downLevel_, float maxDist_, 
									float minDist_,
									int dim_, int cO_, int cB_, boolean tS_, 
                                    String sk_, Vector3f[] sd_, boolean [][][] pm_,
									ViewUserInterface UI_, ViewJProgressBar bar_) {
		image = image_;
		mask = mask_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		upLevel = upLevel_;
        downLevel = downLevel_;
		maxDistance = maxDist_;
        dimension = dim_;
        // mirror connectivity: no!
        connectBack = cB_;
        connectObj = cO_;
		minDistance = minDist_;
		
		thresholdStop = tS_;
        inputSkel = sk_;
        inputSeed = sd_;
        inputPaint = pm_;
        
        UI = UI_;
        progressBar = bar_;
		if(UI_==null||bar_==null){verbose=false;debug=false;}
		// init all the arrays
		try {
			dist = new float[nx][ny][nz];
			geom = new float[nx][ny][nz];
			label = new byte[nx][ny][nz];
			critical = new byte[nx][ny][nz];
			boundary = new BinaryTree(nx*ny*nz, 3, BinaryTree.MINTREE);			
		} catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;

		// init values
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					geom[x][y][z] = 0.0f;
					label[x][y][z] = BACKGROUND;						
					critical[x][y][z] = 0;						
					dist[x][y][z] = 0.0f;
				}
			}
		}
        
		if (dimension==2) dz = 0;
		else if (dimension==3) dz = 1;
		
		if ( (connectObj==4) || (connectObj==6) ) obj = 2;
		else if ( (connectObj==8) || (connectObj==18) ) obj = 3;
		else if (connectObj==26) obj = 4;
		
		if ( (connectBack==4) || (connectBack==6) ) back = 2;
		else if ( (connectBack==8) || (connectBack==18) ) back = 3;
		else if (connectBack==26) back = 4;

		if (debug) MedicUtilPublic.displayMessage("TP:initialisation\n");
	}

	public void finalize() {
		image = null;
		dist = null;
		label = null;
		critical = null;
		geom = null;
		boundary = null;
		System.gc();
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
		label = null;
		//geom = null;
		boundary = null;
		System.gc();
	}
		
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
    public final void setUpLevel(float lv) { upLevel = lv; }
	public final void setDownLevel(float lv) { downLevel = lv; }
	public final void setSkeletonType(String sk) { inputSkel = sk; }
	public final void setImage(float [][][] img) { image = img; }
	
	public final void copyDistanceToImage() { 
		int x,y,z;
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++)
			image[x][y][z] = dist[x][y][z]; 
	}
	
	/**
	 *  higher neighborhood constraint
	 */
	private final float higherNeighborConstraint(int x, int y, int z) {
		float higher;
		int     i,j,l,xi,yj,zl;

		higher = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] == OBJECT) && (dist[xi][yj][zl] > higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  higher neighborhood constraint
	 */
	private final float higherNeighborConstraint(float val, int x, int y, int z) {
		float higher;
		int     i,j,l,xi,yj,zl;

		higher = val;
		//for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < 4) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] == OBJECT) && (dist[xi][yj][zl] > higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  higher neighborhood constraint
	 */
	private final float higherNeighborBidirConstraint(float val, int x, int y, int z) {
		float higher;
		int     i,j,l,xi,yj,zl;

		higher = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] == BGSET) && (dist[xi][yj][zl] > higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborConstraint(int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborBidirConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl]==OBJSET) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBackgroundConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  higher neighborhood constraint
	 */
	private final float higherBoundaryConstraint(int x, int y, int z) {
		float higher,val;
		int     i,j,l,xi,yj,zl;

		val = dist[x][y][z];
		higher = INF;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( ( (label[xi][yj][zl] == BOUNDARY) || (label[xi][yj][zl] == CRITICAL) )
				&& (dist[xi][yj][zl] > val) && (dist[xi][yj][zl] < higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBoundaryConstraint(int x, int y, int z) {
		float lower,val;
		int     i,j,l,xi,yj,zl;

		lower = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
			if (i*i+j*j+l*l < obj) {
				xi=x+i; yj=y+j; zl=z+l; 
				//if ( (label[xi][yj][zl] == BOUNDARY) && (dist[xi][yj][zl] > lower) )
				//	lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == CRITICAL) && (dist[xi][yj][zl] < lower) )
					lower = dist[xi][yj][zl];
			}
		}
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBoundaryNeighbor(int x, int y, int z) {
		float lower,val;
		int     i,j,l,xi,yj,zl;

		lower = -INF;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
			if ( (i*i+j*j+l*l < 4) && (i*i+j*j+l*l > 0) ){
				xi=x+i; yj=y+j; zl=z+l; 
				if ( (label[xi][yj][zl] == BOUNDARY) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == CRITICAL) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == BACKGROUND) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
			}
		}
		if (lower == -INF) lower = dist[x][y][z]-0.01f;
		
		return lower;
	}

	/**
	 *  bigger neighborhood constraint
	 */
	private final int neighborCount(int x, int y, int z) {
		int   val;
		int     i,j,l,xi,yj,zl;

		val = 0;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if (label[xi][yj][zl] >= OBJECT)
				val ++;
		}
	   
		return val;
	}
	private final float neighborFactor(int x, int y, int z) {
		float   val;
		int     i,j,l,xi,yj,zl;

		val = 0.0f;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if (label[xi][yj][zl] >= OBJECT)
				val++;
		}
	   
		return val/(float)connectObj;
	}

	/**
	 *  bigger neighborhood constraint
	 */
	private final float smoothedImage(float[][][] img, int x, int y, int z, float amount) {
		float   val,den;
		int     i,j,l,xi,yj,zl;

		val = img[x][y][z];
		den = 1.0f;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if (label[xi][yj][zl] >= OBJECT) {
				val += amount*img[xi][yj][zl];
				den += amount;
			}
		}
	   
		return val/den;
	}

	/**
	 *  bigger neighborhood constraint
	 */
	private final float closestFactor(int x, int y, int z) {
		float   val;
		int     i,j,l,xi,yj,zl;

		val = 0.0f;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (geom[xi][yj][zl]>val) )
				val = geom[xi][yj][zl];
				
		}
	   
		return val;
	}

	final private float curvatureFactor(int x, int y, int z) {
		float   	num, den, curv;
		float   	ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
	
		// compute the approximate derivatives
		/*
		ux = 0.5f*( dist[x+1][y][z] -dist[x-1][y][z] );
		uy = 0.5f*( dist[x][y+1][z] -dist[x][y-1][z] );
		uz = 0.5f*( dist[x][y][z+1] -dist[x][y][z-1] );
		
		uxx = ( dist[x+1][y][z] -2.0f*dist[x][y][z] +dist[x-1][y][z] );
		uyy = ( dist[x][y+1][z] -2.0f*dist[x][y][z] +dist[x][y-1][z] );
		uzz = ( dist[x][y][z+1] -2.0f*dist[x][y][z] +dist[x][y][z-1] );
		uxy = 0.25f*( dist[x+1][y+1][z] + dist[x-1][y-1][z] 
				    - dist[x-1][y+1][z] - dist[x+1][y-1][z] );
		uyz = 0.25f*( dist[x][y+1][z+1] + dist[x][y-1][z-1] 
				    - dist[x][y-1][z+1] - dist[x][y+1][z-1] );
		uzx = 0.25f*( dist[x+1][y][z+1] + dist[x-1][y][z-1] 
				    - dist[x+1][y][z-1] - dist[x-1][y][z+1] );
		*/
		
		ux = 0.5f*( image[x+1][y][z] -image[x-1][y][z] );
		uy = 0.5f*( image[x][y+1][z] -image[x][y-1][z] );
		uz = 0.5f*( image[x][y][z+1] -image[x][y][z-1] );
		
		uxx = ( image[x+1][y][z] -2.0f*image[x][y][z] +image[x-1][y][z] );
		uyy = ( image[x][y+1][z] -2.0f*image[x][y][z] +image[x][y-1][z] );
		uzz = ( image[x][y][z+1] -2.0f*image[x][y][z] +image[x][y][z-1] );
		uxy = 0.25f*( image[x+1][y+1][z] + image[x-1][y-1][z] 
				    - image[x-1][y+1][z] - image[x+1][y-1][z] );
		uyz = 0.25f*( image[x][y+1][z+1] + image[x][y-1][z-1] 
				    - image[x][y-1][z+1] - image[x][y+1][z-1] );
		uzx = 0.25f*( image[x+1][y][z+1] + image[x-1][y][z-1] 
				    - image[x+1][y][z-1] - image[x-1][y][z+1] );

		num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
		den = (float)Math.sqrt(ux*ux + uy*uy + uz*uz);
		den = 2.0f*den*den*den;
		
		if (den>0.0f) {
			curv = num/den;
		} else {
			curv = 0.0f;
		}
	
		return curv;
	}

	/**
	 * recursive visit of every critical point linked to the current one
	 */
	/*
	private final void visitCritical(int x, int y, int z, float val) {
		int[] vec = new int[3];
		int xi,yj,zl;
		
		label[x][y][z] = BOUNDARY;
		dist[x][y][z] = val;
		vec[0] = x; vec[1] = y; vec[2] = z;
		boundary.addValue(dist[x][y][z], vec);
			
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
			xi=x+i; yj=y+j; zl=z+l;
			if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
				//visitCritical(xi,yj,zl,val);
				label[xi][yj][zl] = BOUNDARY;
				dist[xi][yj][zl] = val;
				vec[0] = xi; vec[1] = yj; vec[2] = zl;
				boundary.addValue(dist[xi][yj][zl], vec);
			
			}
		}
		return;
	}
	 **/

	/**
	 * recursive visit of every critical point linked to the current one
	 */
	/*
	private final void visitCritical(int x, int y, int z, float val) {
		int xi,yj,zl;
		
		dist[x][y][z] = Math.min(val,dist[x][y][z]);
		label[x][y][z] = OBJECT;
			
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
			xi=x+i; yj=y+j; zl=z+l;
			if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
				visitCritical(xi,yj,zl,dist[x][y][z]);
			}
		}
		return;
	}
	*/
	/**
	 * recursive visit of every critical point linked to the current one
	 */
	private final float lowestCritical(int x, int y, int z, float val) {
		int[] vec = new int[3];
		int xi,yj,zl;
		
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
			xi=x+i; yj=y+j; zl=z+l;
			if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
				if (dist[xi][yj][zl] < val) {
					val = lowestCritical(xi,yj,zl,dist[xi][yj][zl]);
				}
			}
		}
		return val;
	}

	/**
	 *  regularity check
	 */
	/*
	private final boolean isRegularObjectPoint(byte[][][] img, int x, int y, int z) {
		boolean regular;
		int Np=0,Nn=0;
	
		//if (debug) MedicUtilPublic.displayMessage("*");
		// counting: check for positive or equal (inside)
		// vs. strictly negative (outside)
		Np = NeighborConnectivity.countNeighborhoods(img,x,y,z,NeighborConnectivity.SUPEQUAL,obj,dz);
		Nn = NeighborConnectivity.countNeighborhoods(img,x,y,z,NeighborConnectivity.INFERIOR,back,dz);
		
		// track possible junction points and local maxima
        if (Np > 1) regular = false;
        else if (Nn > 1) regular = false;
        // add something to avoid 1-pixel holes ??
        else if (Np==0) regular = false;
        else if (Nn==0) regular = false;
        else regular = true;
	
	    return regular;
    }
	private final boolean isRegularBackgroundPoint(byte[][][] img, int x, int y, int z) {
		boolean regular;
		int Np=0,Nn=0;
	
		//if (debug) MedicUtilPublic.displayMessage(":");
		// counting: check for positive or equal (inside)
		// vs. strictly negative (outside)
		Np = NeighborConnectivity.countNeighborhoods(img,x,y,z,NeighborConnectivity.SUPERIOR,back,dz);
		Nn = NeighborConnectivity.countNeighborhoods(img,x,y,z,NeighborConnectivity.INFEQUAL,obj,dz);

		// track possible junction points and local maxima
        if (Np > 1) regular = false;
        else if (Nn > 1) regular = false;
        // add something to avoid 1-pixel holes ??
        else if (Np==0) regular = false;
        else if (Nn==0) regular = false;
        else regular = true;
	
	    return regular;
    }
    */
	/**
	 *  regularity check
	 */
	
	private final boolean isRegularObjectPoint(float[][][] img, int x, int y, int z) {
		byte[][][] template = new byte[3][3][3];
        boolean isRegular;
        
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
            if (img[x+i][y+j][z+l] > img[x][y][z]) template[i+1][j+1][l+1] = 1;
            else if (img[x+i][y+j][z+l] < img[x][y][z]) template[i+1][j+1][l+1] = -1;
            else template[i+1][j+1][l+1] = 0;
        }
        isRegular = isRegularObjectPoint(template,1,1,1);
        template = null;
        return isRegular;
    }
 
	public final boolean isRegularObjectPoint(byte[][][] img, int x, int y, int z) {
		boolean regular;
		int Np=0,Nn=0;
	
		// counting: check for positive or equal (inside)
		// vs. strictly negative (outside)
		if ( (connectObj==6) & (connectBack==26) ) {
			Np = NeighborConnectivity.positiveEqual626Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative26Neighborhood3D(img, x, y, z);
			if (Nn!=1) return false;
        } else if ( (connectObj==6) & (connectBack==18) ) {
			Np = NeighborConnectivity.positiveEqual618Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative18Neighborhood3D(img, x, y, z);
			if (Nn!=1) return false;
        } else if ( (connectObj==6) & (connectBack==6) ) {
			Np = NeighborConnectivity.positive626Neighborhood3D(img, x, y, z);
            if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual26Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
			Np = NeighborConnectivity.positive26Neighborhood3D(img, x, y, z);
            if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual626Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        } else if (connectObj==18) {
			Np = NeighborConnectivity.positiveEqual18Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative618Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        } else if (connectObj==26) {
			Np = NeighborConnectivity.positiveEqual26Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative626Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        } else if (connectObj==4) {
			Np = NeighborConnectivity.positiveEqual4Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative8Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        } else if (connectObj==8) {
			Np = NeighborConnectivity.positiveEqual8Neighborhood3D(img, x, y,z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negative4Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        }
		return true;
    }
	public final boolean isRegularBackgroundPoint(byte[][][] img, int x, int y, int z) {
		boolean regular;
		int Np=0,Nn=0;
	
		// counting: check for strictly positive(inside)
		// vs. negative or equal (outside)
		if ( (connectBack==6) & (connectObj==26) ) {
			Np = NeighborConnectivity.positive626Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual26Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        } else if ( (connectBack==6) & (connectObj==18) ) {
			Np = NeighborConnectivity.positive618Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual18Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        } else if ( (connectBack==6) & (connectObj==6) ) {
			Np = NeighborConnectivity.positive626Neighborhood3D(img, x, y, z);
            if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual26Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
			Np = NeighborConnectivity.positive26Neighborhood3D(img, x, y, z);
            if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual626Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;            
        } else if (connectBack==18) {
			Np = NeighborConnectivity.positive18Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual618Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        } else if (connectBack==26) {
			Np = NeighborConnectivity.positive26Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual626Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        } else if (connectBack==4) {
			Np = NeighborConnectivity.positive4Neighborhood3D(img, x, y, z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual8Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        } else if (connectBack==8) {
			Np = NeighborConnectivity.positive8Neighborhood3D(img, x, y,z);
			if (Np!=1) return false;
            Nn = NeighborConnectivity.negativeEqual4Neighborhood3D(img, x, y, z);
            if (Nn!=1) return false;
        }
        
        return true;
    }
	
	/**
	 *  initialize the labels
	 */
	public final void initUpwardGeometricLabels() {
        int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( image[x][y][z] <= upLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (image[x+i][y+j][z+l] <= upLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardGeometricLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+downLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( image[x][y][z] >= downLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (image[x+i][y+j][z+l] >= downLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = -1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
    
	public final void initUpwardSmoothingLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( image[x][y][z] <= upLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
						geom[x][y][z] = 1.0f;
						dist[x][y][z] = image[x][y][z];
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (image[x+i][y+j][z+l] <= upLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initUpwardSeedLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from seed points and up
        for (int n=0;n<inputSeed.length;n++) {
            // initial object 
            int x = (int)inputSeed[n].X;
            int y = (int)inputSeed[n].Y;
            int z = (int)inputSeed[n].Z;
            label[x][y][z] = OBJECT;
            geom[x][y][z] = 1.0f;
            dist[x][y][z] = image[x][y][z];
        }
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if (label[x][y][z]==OBJECT) {
                        // do nothing
                    } else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (label[x+i][y+j][z+l] == OBJECT) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initUpwardPaintLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( inputPaint[x][y][z] ) {
						// initial object 
						label[x][y][z] = OBJECT;
						geom[x][y][z] = 1.0f;
						dist[x][y][z] = image[x][y][z];
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (inputPaint[x+i][y+j][z+l]) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
    
	public final void initDownwardSmoothingLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( image[x][y][z] >= downLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 1.0f;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (image[x+i][y+j][z+l] >= downLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
							//boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardSeedLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from seed points and up
        for (int n=0;n<inputSeed.length;n++) {
            // initial object 
            int x = (int)inputSeed[n].X;
            int y = (int)inputSeed[n].Y;
            int z = (int)inputSeed[n].Z;
            label[x][y][z] = OBJECT;
            geom[x][y][z] = 1.0f;
            dist[x][y][z] = image[x][y][z];
        }
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( label[x][y][z] == OBJECT ) {
						// do nothing
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (label[x+i][y+j][z+l] == OBJECT) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            //boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardPaintLabels() {
		int[] vec = new int[3];
		if (debug) MedicUtilPublic.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( inputPaint[x][y][z] ) {
						// initial object 
						label[x][y][z] = OBJECT;
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 1.0f;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (inputPaint[x+i][y+j][z+l]) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            //boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
			
	/**
	 *  propagate distances from a contour to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		initUpwardGeometricLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax, false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			higher = higherNeighborConstraint(x, y, z);
            
			if (higher>val) {
				val = higher;
                dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);
				critical[x][y][z] += 1;
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
               // test for update
                label[x][y][z] = TEST;

                if (isRegularBackgroundPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val+1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
               }
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val+1.0f;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateUpwardGeometricDistance

	/**
	 *  propagate distances from a contour to the outside (inferior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		initDownwardGeometricLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(x, y, z);
            
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
                 vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);               
				critical[x][y][z] -= 1;
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val-1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
               }
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val-1.0f;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardGeometricDistance

	/**
	 *  propagate distances from a contour to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardUnconstrainedGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		initUpwardGeometricLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			higher = higherNeighborConstraint(x, y, z);
            
			if (higher>val) {
				val = higher;
                dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);
				critical[x][y][z] += 1;
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
               // test for update
                label[x][y][z] = TEST;

                //if (isRegularBackgroundPoint(label,x,y,z) ) {
                if ( true ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val+1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
               }
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val+1.0f;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateUpwardGeometricDistance

	/**
	 *  propagate distances from a contour to the outside (inferior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardUnconstrainedGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		initDownwardGeometricLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(x, y, z);
            
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
                 vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);               
				critical[x][y][z] -= 1;
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                //if (isRegularObjectPoint(label,x,y,z) ) {
                if ( true ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val-1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
               }
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val-1.0f;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardGeometricDistance

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardExactSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			
			boundary.removeFirst();

			//val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if (isRegularBackgroundPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                        //if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.max(image[xi][yj][zl], val+minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						//if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val+minDistance, vec);
                            //boundary.addValue(val, vec);
							//followCriticalPath(xi,yj,zl,val+minDistance);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardExactSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardExactSmoothing() {
	
		int 		t, stop;
		float   	val, dival, prec, rank;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
		else initDownwardSmoothingLabels();
		if (Boolean.getBoolean("macruise.topology.trace")) {
			System.err.println("TRACE_INIT dimensions="+nx+","+ny+","+nz+" obj="+obj+" back="+back+" seed=5,5,5");
			for (int xx=4;xx<=6;xx++) for (int yy=4;yy<=6;yy++) for (int zz=4;zz<=6;zz++)
				System.err.println("TRACE_LABEL "+xx+","+yy+","+zz+"="+label[xx][yy][zz]+" dist="+dist[xx][yy][zz]);
			if (boundary.isNotEmpty())
				System.err.println("TRACE_QUEUE_HEAD "+boundary.getFirstIndex(0)+","+boundary.getFirstIndex(1)+","+boundary.getFirstIndex(2)+" value="+boundary.getFirst());
		}
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		dival = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			boolean traceThisEvent = Boolean.getBoolean("macruise.topology.trace");
			int traceEventId = diagnosticQueueEvents++;
			if (traceThisEvent)
				System.err.println("TRACE_EVENT "+traceEventId+" point="+x+","+y+","+z+" value="+val+" lower="+lower+" label="+label[x][y][z]);
			if (lower<val) {
				val = lower;
				//dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				boolean regular = isRegularObjectPoint(label,x,y,z);
				if (traceThisEvent) {
					int foreground = NeighborConnectivity.positiveEqual618Neighborhood3D(label,x,y,z);
					int background = NeighborConnectivity.negative18Neighborhood3D(label,x,y,z);
					System.err.println("TRACE_TEST "+traceEventId+" foreground="+foreground+" background="+background+" regular="+regular);
				}

                if (regular) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            //boundary.addValue(val-minDistance, vec);
                            boundary.addValue(val, vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = Math.min(image[xi][yj][zl], val-minDistance);
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(Math.min(image[xi][yj][zl], val-minDistance), vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}				
				if (traceThisEvent) System.err.println("TRACE_DECISION "+traceEventId+" status="+label[x][y][z]+" distance="+dist[x][y][z]);
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		if (Boolean.getBoolean("macruise.topology.trace"))
			System.err.println("TRACE_EVENT_TOTAL "+diagnosticQueueEvents);
		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z] == CRITICAL) ) {
				if (thresholdStop) dist[x][y][z] = upLevel;
				else dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateDownwardGeometryDrivenSmoothing() {
	
		int 		t, stop;
		float   	val, dival, prec, rank;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		// init the labels from the threshold and up
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( image[x][y][z] >= downLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 0.0f;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (image[x+i][y+j][z+l] >= downLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(geom[x][y][z], vec);
                            //boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		dival = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			rank = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
				dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(rank+1, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
                dist[x][y][z] = val;
				
                if (isRegularObjectPoint(dist,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(rank+1, vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(rank+1, vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
                    dist[x][y][z] = 0.0f;
				}				
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z] == CRITICAL) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing
	*/
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateBidirectionalExactSmoothing() {
	
		int 		t, stop;
		float   	val, prec, rank;
		float   	lower,higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		boolean		isProcessed, isRegular;
		float		upMax,downMin;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		// init the labels from the threshold and up
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
				label[x][y][z] = BGSET;
			} else if ( mask[x][y][z]==false ) {
				label[x][y][z] = BGSET;
			} else if ( image[x][y][z] <= upLevel ) {
				// initial background 
				label[x][y][z] = BGSET;
				dist[x][y][z] = image[x][y][z];
			} else if ( image[x][y][z] >= downLevel ) {
				// initial object 
				label[x][y][z] = OBJSET;
				dist[x][y][z] = image[x][y][z];
			} else {
				// boundary or background ?
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
					if (i*i+j*j+l*l < obj) {
						if ( (image[x+i][y+j][z+l] >= downLevel) && (mask[x+i][y+j][z+l]) ) {
							label[x][y][z] = OBJBOUND;
						}					
					}
					if (i*i+j*j+l*l < back) {
						if ( (image[x+i][y+j][z+l] <= upLevel) && (mask[x+i][y+j][z+l]) ) {
							label[x][y][z] = BGBOUND;
						}
						if (!mask[x+i][y+j][z+l]) 
							label[x][y][z] = BGBOUND;
						if ( (x+i==0) || (y+j==0) || (z+l==0) || (x+i==nx-1) || (y+j==ny-1) || (z+l==nz-1) ) {
							label[x][y][z] = BGBOUND;
						}
					}
				}
				if (label[x][y][z] == OBJBOUND) {
					// object boundary: stack to the tree
					dist[x][y][z] = image[x][y][z];
					vec[0] = x; vec[1] = y; vec[2] = z;
                    boundary.addValue(downLevel-dist[x][y][z], vec);
				} else if (label[x][y][z] == BGBOUND) {
					// background boundary: stack to the tree
					dist[x][y][z] = image[x][y][z];
					vec[0] = x; vec[1] = y; vec[2] = z;
                    boundary.addValue(dist[x][y][z]-upLevel, vec);
				} else {
					label[x][y][z] = UNSET;
					dist[x][y][z] = 0.0f;
				}
			}
		}
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		downMin = downLevel;
		upMax = upLevel;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			isProcessed = false;
			if (label[x][y][z]==OBJBOUND) {
				val = downLevel - val;
                if (val < downMin) downMin = val;
				lower = lowerNeighborBidirConstraint(val, x, y, z);
				if (lower<val) {
					// set as critical ?
					label[x][y][z] = OBJCRIT;
					isProcessed = false;
				} else {
					isProcessed = true;
				}
			} else if (label[x][y][z]==BGBOUND) {
				val = upLevel + val;
                if (val > upMax) upMax = val;
				higher = higherNeighborBidirConstraint(val, x, y, z);
				if (higher<val) {
					label[x][y][z] = BGCRIT;
					isProcessed = false;
				} else {
					isProcessed = true;					
				}
			}
            if (isProcessed) {
				isRegular = false;
				if (label[x][y][z]==OBJBOUND) {
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BGSET) && (dist[xi][yj][zl]>=val) ) {
                            label[xi][yj][zl] = OBJBOUND;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(downLevel-val, vec);
                        }
                    }
					// test for update
					label[x][y][z] = OBJTEST;
					if (isRegularObjectPoint(label,x,y,z) ) isRegular = true;
					else isRegular = false;
				} else if (label[x][y][z]==BGBOUND) {
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == OBJSET) && (dist[xi][yj][zl]<=dist[x][y][z]) ) {
                            label[xi][yj][zl] = BGBOUND;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val-upLevel, vec);
                        }
                    }
					// test for update
					label[x][y][z] = BGTEST;
					if (isRegularBackgroundPoint(label,x,y,z) ) isRegular = true;
					else isRegular = false;
				} 

				if (isRegular) {
                    // regular point
					if (label[x][y][z]==OBJTEST) {
						label[x][y][z] = OBJSET;
						dist[x][y][z] = val;

						// search for new or critical neighbors
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							xi=x+i; yj=y+j; zl=z+l;
							if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == OBJCRIT) ) {
								label[xi][yj][zl] = OBJBOUND;
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(downLevel-val, vec);
							}
							if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == UNSET) ) {
								label[xi][yj][zl] = OBJBOUND;
								dist[xi][yj][zl] = image[xi][yj][zl];
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(downLevel-image[xi][yj][zl], vec);
							}
                            if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BGCRIT) ) {
								label[xi][yj][zl] = OBJBOUND;
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(downLevel-val, vec);
							}
						}
					} else if (label[x][y][z]==BGTEST) {
						label[x][y][z] = BGSET;
						dist[x][y][z] = val;

						// search for new or critical neighbors
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							xi=x+i; yj=y+j; zl=z+l;
							if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BGCRIT) ) {
								label[xi][yj][zl] = BGBOUND;
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(val-upLevel, vec);
							}
							if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == UNSET) ) {
								label[xi][yj][zl] = BGBOUND;
								dist[xi][yj][zl] = image[xi][yj][zl];
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(image[xi][yj][zl]-upLevel, vec);
							}
                            if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == OBJCRIT) ) {
								label[xi][yj][zl] = BGBOUND;
								vec[0] = xi; vec[1] = yj; vec[2] = zl;
								boundary.addValue(val-upLevel, vec);
							}
						}
					}
				} else {
					// critical point
					if (label[x][y][z]==OBJTEST) {
						label[x][y][z] = OBJCRIT;
					} else if (label[x][y][z]==BGTEST) {
						label[x][y][z] = BGCRIT;
					}
				}				
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if (label[x][y][z] == BGCRIT) {
				dist[x][y][z] = downLevel;
			} else if (label[x][y][z] == OBJCRIT) {
				dist[x][y][z] = upLevel;
			} else if (label[x][y][z]==UNSET) {
				dist[x][y][z] = 0.5f*(downLevel+upLevel);
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing
	*/
	
    public final void updateCriticalTree(int x, int y, int z, float val) {
        // search for all critical neighbors: if none, put the point into the boundary
        boolean hasNeighbors = false;
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
            if ( (i*i+j*j+l*l < back) && (label[x+i][y+j][z+l] == CRITICAL) ) {
                hasNeighbors = true;
                label[x][y][z] = MASK;
                updateCriticalTree(x+i,y+j,z+l,val);
             }
        }
        hasNeighbors = false;
        if (!hasNeighbors) {    
            label[x][y][z] = BOUNDARY;
            //dist[xi][yj][zl] = val;
            int[] vec = new int[3]; 
            vec[0] = x; vec[1] = y; vec[2] = z;
            boundary.addValue(val, vec);
            vec = null;
        } else {
            // change the label back to critical
            label[x][y][z] = CRITICAL;
        }
    }
    
    
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateDownwardCriticalSmoothing() {
	
		int 		t, stop;
		float   	val, dival, prec, rank;
		float   	lower,higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		dival = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
                //dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            /*
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val, vec);
                             //
                            updateCriticalTree(xi,yj,zl,val);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(image[xi][yj][zl], vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}				
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");
        /*
        // init the labels from the threshold and up
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( label[x][y][z] == OBJECT ) {
						// initial object : reset
						label[x][y][z] = OBJECT;
						geom[x][y][z] = 1.0f;
						dist[x][y][z] = image[x][y][z];
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (image[x+i][y+j][z+l] <= upLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}

/*
        for (int N=0;N<0;N++) {
        // try to propagate starting on the background points not reached
        boundary.reset();
		boundary.setMaxTree();
        // init the labels from the threshold and up
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( label[x][y][z] == OBJECT ) {
						// initial object 
					} else if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z]==CRITICAL) ) {			
						label[x][y][z] = BOUNDARY;
						// boundary: stack to the tree
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 1.0f;
						vec[0] = x; vec[1] = y; vec[2] = z;
                        //boundary.addValue(geom[x][y][z], vec);
                        boundary.addValue(dist[x][y][z], vec);                        
					} else if (label[x][y][z] == CRITICAL) {
                        label[x][y][z] = BACKGROUND;
                        dist[x][y][z] = 0.0f;
                    }
				}
			}
		}

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		dival = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
                //dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val, vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(image[xi][yj][zl], vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}				
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");
        }
        
		// make a new tree going up from the boundaries
		// re-init labels
		boundary.reset();
		boundary.setMinTree();
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
				label[x][y][z] = FRONTIER;
			} else if ( mask[x][y][z]==false ) {
				label[x][y][z] = MASK;
			} else if ( (label[x][y][z] == CRITICAL) || (label[x][y][z]==BACKGROUND) ) {
				// critical object : check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
					if (i*i+j*j+l*l < 4) {
						if ( (label[x+i][y+j][z+l]==FRONTIER) || (label[x+i][y+j][z+l]==MASK) ) {
								label[x][y][z] = MBOUND;
						}					
					}
				}
				if (label[x][y][z] == MBOUND) {
					// boundary: stack to the tree
					dist[x][y][z] = lowerNeighborConstraint(image[x][y][z],x,y,z);
					geom[x][y][z] = 1.0f;
					vec[0] = x; vec[1] = y; vec[2] = z;
                    boundary.addValue(dist[x][y][z], vec);
				} else {
					label[x][y][z] = MBACK;
				}
			} else if (label[x][y][z]==BACKGROUND) {
                label[x][y][z]=MBOUND;
                dist[x][y][z] = image[x][y][z];
				geom[x][y][z] = 1.0f;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);
			} else if (label[x][y][z]==OBJECT) {
                label[x][y][z]=OBJECT;
			} 
		}
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
            
			if (lower<val) {
				val = lower;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else {
				// test for update
                label[x][y][z] = MTEST;

                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = MOBJ;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == MBACK) ) {
                            label[xi][yj][zl] = MBOUND;
                            dist[xi][yj][zl] = image[xi][yj][zl];
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = MBOUND;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val, vec);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
//
		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				//dist[x][y][z] = 0.0f;
                dist[x][y][z] = lowerNeighborConstraint(image[x][y][z],x,y,z);
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing
	*/
	
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateDownwardLoopedSmoothing() {
	
		int 		t, stop;
		float   	val, dival, prec, rank;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();
			
			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(image[xi][yj][zl], vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
					// replace value with closest
					val  = lowerBoundaryNeighbor(x,y,z);
					//val = dist[x][y][z] - 0.001f;
					dist[x][y][z] = val;
                    vec[0] = x; vec[1] = y; vec[2] = z;
                    boundary.addValue(val, vec);
                }
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z] == CRITICAL) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing
	*/
	
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardApproxSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			boundary.removeFirst();

			//val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if ( (val >= image[x][y][z] + maxDistance) || (isRegularBackgroundPoint(label,x,y,z) ) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = Math.min(val, image[x][y][z] + maxDistance);
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.max(image[xi][yj][zl],val+minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						// variant
						if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val+minDistance, vec);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = Math.min(image[x][y][z]+maxDistance,val);
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardApproxSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardApproxSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax,false);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
            if (lower<val) {
				// bound if necessary
				val = lower;
                //dist[x][y][z] = val;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                label[x][y][z] = TEST;
				
                if ( (val <= image[x][y][z] - maxDistance) || (isRegularObjectPoint(label,x,y,z) ) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = Math.max(val, image[x][y][z] - maxDistance);
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.min(image[xi][yj][zl],val-minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						// variant: keep the critical inside with lower value
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val-minDistance, vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if (label[x][y][z] == CRITICAL) {
				dist[x][y][z] = Math.max(image[x][y][z]-maxDistance, val);
			} else if (label[x][y][z] == BACKGROUND) {
				dist[x][y][z] = image[x][y][z];
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardApproxSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateUpwardBoostedSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			
			boundary.removeFirst();

			val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(geom[x][y][z]*val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if (isRegularBackgroundPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(geom[x][y][z]*val, vec);
                        }
						if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							geom[xi][yj][zl] = (1.0f+maxDistance)*geom[x][y][z];
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(geom[xi][yj][zl]*dist[xi][yj][zl], vec);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardExactSmoothing
	*/
	
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateDownwardBoostedSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[]		vec = new int[3];
		float[]		neighbors = new float[26];
		int			nNeigh;
		float			maxNeigh;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
            
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(geom[x][y][z]*dist[x][y][z], vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                if (isRegularObjectPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(geom[xi][yj][zl]*dist[xi][yj][zl], vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
							geom[xi][yj][zl] = (1.0f-maxDistance)*geom[x][y][z];
							dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
							boundary.addValue(geom[xi][yj][zl]*dist[xi][yj][zl], vec);
						}
					}
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val;
			}
		}
		// set the critical points : the heavy technique
		// note: to be optimal, we should make a tree with all endpoints
		// ranked by decreasing value -> propagate all at once
		int Ncrit = 0, Nprec = 0;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				Ncrit++;
			}
		}
		while ( (Ncrit>0) && (Ncrit!=Nprec) ) {		
			if (debug) MedicUtilPublic.displayMessage("critical: "+Ncrit+"\n");
			Nprec = Ncrit;
			for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
				if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
					label[x][y][z] = TEST;
					if (isRegularObjectPoint(label,x,y,z)) {
						dist[x][y][z] = lowerBackgroundConstraint(image[x][y][z],x,y,z);
						label[x][y][z] = OBJECT;
						Ncrit--;
						critical[x][y][z] = 2;
					} else {
						label[x][y][z] = CRITICAL;
					}
				}
			}
		}
		
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardBoostedSmoothing
	*/
	
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateUpwardJoinedSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			
			boundary.removeFirst();

			//val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                //if (isRegularBackgroundPoint(label,x,y,z) ) {
                if (true) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val, vec);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardJoinedSmoothing
	 **/

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	/*
	public final void propagateDownwardJoinedSmoothing() {
	
		int 		t, stop;
		float   	val, prec, rank;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) MedicUtilPublic.displayMessage("TP:start upward distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();
		
		if (debug) MedicUtilPublic.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (verbose) && (t%mod==0) ) {
				progressBar.updateValue(t*100/tmax);
			}
			if ( (debug) && (t%10000==0)) MedicUtilPublic.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
                //dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                //if (isRegularObjectPoint(label,x,y,z) ) {
                if (true) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val, vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = image[xi][yj][zl];
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(image[xi][yj][zl], vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}				
            }
		}
		if (debug) MedicUtilPublic.displayMessage("TP:critical points\n");

		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z] == CRITICAL) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) MedicUtilPublic.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardJoinedSmoothing
	 */

    /** 
	 *	returns some internal functions 
	 */
	public final float[][][] exportDistance() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = dist[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances
	public final float[][][] exportGeometry() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = geom[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances
	public final float[][][] exportLabel() {
		int 	x,y,z,k;
		float[][][]	lbs = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					lbs[x][y][z] = (float)label[x][y][z];
				}
			}
		}
		return lbs;
	} // exportDistances
	public final float[][][] exportCriticalPoints() {
		int 	x,y,z,k;
		float[][][]	cpt = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					cpt[x][y][z] = (float)critical[x][y][z];
				}
			}
		}
		return cpt;
	} // exportDistances
	public final float[][][] exportImage() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = image[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances

	public final void visitCriticalPoint(int x, int y, int z,float val) {
		
		// endpoint
		label[x][y][z] = OBJECT;
		//dist[x][y][z] = Math.min(lowerNeighborConstraint(x,y,z),val);
		dist[x][y][z] = 0.0f;
		// propagate
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
			if ( (label[x+i][y+j][z+l]==CRITICAL) || (label[x+i][y+j][z+l]==BACKGROUND) ) {
				label[x+i][y+j][z+l] = TEST;
				if (isRegularObjectPoint(label,x+i,y+j,z+l)) {
					visitCriticalPoint(x+i,y+j,z+l,dist[x][y][z]);
				} else {
					label[x+i][y+j][z+l] = CRITICAL;
				}
			}
		}
	} // visitCriticalPoint
}
