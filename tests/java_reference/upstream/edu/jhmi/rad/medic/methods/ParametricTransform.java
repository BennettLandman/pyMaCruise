package edu.jhmi.rad.medic.methods;


import java.io.*;
import java.util.*;
import java.lang.*;

import edu.jhmi.rad.medic.structures.*;
import edu.jhmi.rad.medic.utilities.*;

/**
 *
 *   Centralized generic functions for spatial transformation of coordinates
 *
 *
 *	@version    Feb 2009
 *	@author     Pierre-Louis Bazin
 *		
 *
*/
public class ParametricTransform {
	
	int type;
	private static final	int   	NONE = 0;
	private static final	int   	RIGID = 11;
	private static final	int   	SINGLE_SCALE = 12;
	private static final	int		SCALED_RIGID = 13;
	private static final	int   	QUADRATIC_SCALE = 21;
	private static final	int   	TRI_QUADRATIC_SCALE = 22;
	private static final	int   	FULLY_AFFINE = 31;
	private static final	int   	FULLY_QUADRATIC = 32;
	private static final	int   	FULLY_CUBIC = 33;
	
	float 	x0i, y0i, z0i;	// image center
	float 	rix, riy, riz;	// image resolutions
	int 	nix, niy, niz;	// image dimensions
	
	float 	x0t, y0t, z0t;	// template center
	float 	rtx, rty, rtz;	// template resolutions
	int 	ntx, nty, ntz;	// template dimensions
	
	int 	Nt;				// transform dimensions
	
	float Iscale, Iscale2;
	
	boolean debug = true;
	
	public ParametricTransform(String type_, 
								float x0i_, float y0i_, float z0i_,
								float rix_, float riy_, float riz_,
								int nix_, int niy_, int niz_,
								float x0t_, float y0t_, float z0t_,
								float rtx_, float rty_, float rtz_,
								int ntx_, int nty_, int ntz_) {
		
		if (type_.equals("rigid")) { type = RIGID; Nt = 6; } 
		else if (type_.equals("single_scale")) { type = SINGLE_SCALE; Nt = 7; }
		else if (type_.equals("scaled_rigid")) { type = SCALED_RIGID; Nt = 9; }
		else if (type_.equals("fully_affine")) { type = FULLY_AFFINE; Nt = 12; }
		else if (type_.equals("quadratic_scale")) { type = QUADRATIC_SCALE; Nt = 8; }
		else if (type_.equals("tri_quadratic_scale")) { type = TRI_QUADRATIC_SCALE; Nt = 12; }
		else if (type_.equals("fully_quadratic")) { type = FULLY_QUADRATIC; Nt = 30; }
		else if (type_.equals("fully_cubic")) { type = FULLY_CUBIC; Nt = 60; }
		else { type = NONE; Nt = 0; }
		
		if (debug) MedicUtilPublic.displayMessage("transform type: "+type_+"("+type+")\n");
		
		x0i = x0i_; y0i = y0i_; z0i = z0i_;
		rix = rix_; riy = riy_; riz = riz_;
		nix = nix_; niy = niy_; niz = niz_;
		
		x0t = x0t_; y0t = y0t_; z0t = z0t_;
		rtx = rtx_; rty = rty_; rtz = rtz_;
		ntx = ntx_; nty = nty_; ntz = ntz_;
		
		Iscale2 = 0.25f*( nix*rix*nix*rix + niy*riy*niy*riy + niz*riz*niz*riz );
		Iscale = (float)Math.sqrt(Iscale2);	
	}
	
	public final void finalize() {
	}
	
	public final int getDimension() { return Nt; }
	
	public final String getTransformType() {
		if (type==NONE) 						return "none";
		else if (type==RIGID) 					return "rigid";
		else if (type==SINGLE_SCALE) 			return "single_scale";
		else if (type==SCALED_RIGID) 			return "scaled_rigid";
		else if (type==QUADRATIC_SCALE) 		return "quadratic_scale";
		else if (type==TRI_QUADRATIC_SCALE) 	return "tri_quadratic_scale";
		else if (type==FULLY_AFFINE) 			return "fully_affine";
		else if (type==FULLY_QUADRATIC) 		return "fully_quadratic";
		else if (type==FULLY_CUBIC) 			return "fully_cubic";
		else 									return "not found!";
	}		
	
	/** check if the transform is linear in the image coordinates (and thus pre-computable) */
	public final boolean isLinear() {
		if (type==NONE || type==RIGID || type==SINGLE_SCALE || type==SCALED_RIGID || type==FULLY_AFFINE) return true;
		else return false;
	}
	
	/** check if the transform uses a rotation matrix */
	public final boolean useRotation() {
		if (type==RIGID || type==SINGLE_SCALE || type==SCALED_RIGID || type==QUADRATIC_SCALE || type==TRI_QUADRATIC_SCALE) return true;
		else return false;
	}
	
	/** 
	 *	computes the transformed coordinates from image to template space
	 *	with scaling (2: half res, 4: quarter res, etc)
	 */
	public final float[] imageToTemplate(int x,int y,int z, float[] trans, float[][] rot, float scale) {
		float[] X = new float[3];
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		if (type==NONE) {
			X[0] = ( Xi/rtx + x0t)/scale;
			X[1] = ( Yi/rty + y0t)/scale;
			X[2] = ( Zi/rtz + z0t)/scale;
		} else if (type==RIGID) {
			X[0] = ( (rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( (rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( (rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			
			X[0] = ( factor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( factor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( factor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==SCALED_RIGID) {
			X[0] = ( (1.0f+trans[6]/Iscale)*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( (1.0f+trans[7]/Iscale)*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( (1.0f+trans[8]/Iscale)*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		
		} else if (type==QUADRATIC_SCALE) {
			float norm = Xi*Xi + Yi*Yi + Zi*Zi;
			float factor = 1.0f + trans[6]/Iscale + trans[7]/Iscale*(2.0f*norm/Iscale2-1.0f);
			
			X[0] = ( factor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( factor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( factor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==TRI_QUADRATIC_SCALE) {
			float normx = rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3]; normx *= normx;
			float normy = rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4]; normy *= normy;
			float normz = rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5]; normz *= normz;

			float factorx = 1.0f + trans[6]/Iscale + trans[ 9]/Iscale*(2.0f*normx/Iscale2-1.0f);
			float factory = 1.0f + trans[7]/Iscale + trans[10]/Iscale*(2.0f*normy/Iscale2-1.0f);
			float factorz = 1.0f + trans[8]/Iscale + trans[11]/Iscale*(2.0f*normz/Iscale2-1.0f);
			
			X[0] = ( factorx*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( factory*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( factorz*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==FULLY_AFFINE) {
			X[0] = ( ((1+trans[0])*Xi+   trans[1] *Yi+   trans[2] *Zi+trans[ 9])/rtx + x0t)/scale;
			X[1] = ( (   trans[3] *Xi+(1+trans[4])*Yi+   trans[5] *Zi+trans[10])/rty + y0t)/scale;
			X[2] = ( (   trans[6] *Xi+   trans[7] *Yi+(1+trans[8])*Zi+trans[11])/rtz + z0t)/scale;
		} else if (type==FULLY_QUADRATIC) {
			X[0] = ( ( (1+trans[ 0])*Xi+   trans[ 1] *Yi+   trans[ 2] *Zi + trans[ 9] 
					  +trans[12]*(2.0f*Xi*Xi/Iscale-Iscale)+trans[13]*(2.0f*Yi*Yi/Iscale-Iscale)+trans[14]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[15]*Xi*Yi/Iscale +trans[16]*Yi*Zi/Iscale +trans[17]*Zi*Xi/Iscale )/rtx + x0t)/scale;
			X[1] = ( (    trans[ 3] *Xi+(1+trans[ 4])*Yi+   trans[ 5] *Zi + trans[10] 
					  +trans[18]*(2.0f*Xi*Xi/Iscale-Iscale)+trans[19]*(2.0f*Yi*Yi/Iscale-Iscale)+trans[20]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[21]*Xi*Yi/Iscale +trans[22]*Yi*Zi/Iscale +trans[23]*Zi*Xi/Iscale )/rty + y0t)/scale;
			X[2] = ( (    trans[ 6] *Xi+   trans[ 7] *Yi+(1+trans[ 8])*Zi + trans[11] 
					  +trans[24]*(2.0f*Xi*Xi/Iscale-Iscale)+trans[25]*(2.0f*Yi*Yi/Iscale-Iscale)+trans[26]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[27]*Xi*Yi/Iscale +trans[28]*Yi*Zi/Iscale +trans[29]*Zi*Xi/Iscale )/rtz + z0t)/scale;
		} else if (type==FULLY_CUBIC) {
			X[0] = ( ( (1+trans[ 0])*Xi+   trans[ 1] *Yi+   trans[ 2] *Zi + trans[ 9] 
					  +trans[12]*(2.0f*Xi*Xi/Iscale-Iscale) +trans[13]*(2.0f*Yi*Yi/Iscale-Iscale) +trans[14]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[15]*Xi*Yi/Iscale +trans[16]*Yi*Zi/Iscale +trans[17]*Zi*Xi/Iscale
					  +trans[30]*(4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi) +trans[31]*(4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi) +trans[32]*(4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)
					  +trans[33]*(2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale +trans[34]*(2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale +trans[35]*(2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale 
					  +trans[36]*(2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale +trans[37]*(2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale +trans[38]*(2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale 
					  +trans[39]*Xi*Yi*Zi/Iscale2)/rtx + x0t)/scale;
			X[1] = ( (    trans[ 3] *Xi+(1+trans[ 4])*Yi+   trans[ 5] *Zi + trans[10] 
					  +trans[18]*(2.0f*Xi*Xi/Iscale-Iscale)+trans[19]*(2.0f*Yi*Yi/Iscale-Iscale)+trans[20]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[21]*Xi*Yi/Iscale +trans[22]*Yi*Zi/Iscale +trans[23]*Zi*Xi/Iscale
					  +trans[40]*(4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi) +trans[41]*(4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi) +trans[42]*(4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)
					  +trans[43]*(2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale +trans[44]*(2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale +trans[45]*(2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale 
					  +trans[46]*(2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale +trans[47]*(2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale +trans[48]*(2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale 
					  +trans[49]*Xi*Yi*Zi/Iscale2)/rty + y0t)/scale;
			X[2] = ( (    trans[ 6] *Xi+   trans[ 7] *Yi+(1+trans[ 8])*Zi + trans[11] 
					  +trans[24]*(2.0f*Xi*Xi/Iscale-Iscale)+trans[25]*(2.0f*Yi*Yi/Iscale-Iscale)+trans[26]*(2.0f*Zi*Zi/Iscale-Iscale)
					  +trans[27]*Xi*Yi/Iscale +trans[28]*Yi*Zi/Iscale +trans[29]*Zi*Xi/Iscale
					  +trans[50]*(4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi) +trans[51]*(4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi) +trans[52]*(4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)
					  +trans[53]*(2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale +trans[54]*(2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale +trans[55]*(2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale 
					  +trans[56]*(2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale +trans[57]*(2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale +trans[58]*(2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale 
					  +trans[59]*Xi*Yi*Zi/Iscale2)/rtz + z0t)/scale;
		} 
		return X;
	}

	public final float[][] imageToTemplateDerivatives(int x,int y,int z,float[] trans, float[][] rot, float[][] dRa, float[][] dRb, float[][] dRc, float scale) {
		float[][] dX = new float[3][Nt];
		
		for (int i=0;i<3;i++) 
			for (int t=0;t<Nt;t++)
				dX[i][t] = 0.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
			return null;
		} else if (type==RIGID) {
			dX[0][0] = ( (dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( (dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( (dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = 1.0f/rtx/scale;
			
			dX[1][0] = ( (dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( (dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( (dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = 1.0f/rty/scale;
			
			dX[2][0] = ( (dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( (dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( (dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = 1.0f/rtz/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dX[0][0] = ( factor*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( factor*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( factor*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = factor/rtx/scale;
			dX[0][6] = ( dfactor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][0] = ( factor*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( factor*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( factor*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = factor/rty/scale;
			dX[1][6] = ( dfactor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][0] = ( factor*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( factor*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( factor*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = factor/rtz/scale;
			dX[2][6] = ( dfactor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==SCALED_RIGID) {
			float factorx = 1.0f + trans[6]/Iscale;
			float factory = 1.0f + trans[7]/Iscale;
			float factorz = 1.0f + trans[8]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dX[0][0] = ( factorx*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( factorx*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( factorx*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = factorx/rtx/scale;
			dX[0][6] = ( dfactor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][0] = ( factory*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( factory*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( factory*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = factory/rty/scale;
			dX[1][7] = ( dfactor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][0] = ( factorz*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( factorz*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( factorz*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = factorz/rtz/scale;
			dX[2][8] = ( dfactor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==QUADRATIC_SCALE) {
			float norm = Xi*Xi + Yi*Yi + Zi*Zi;
			float factor = 1.0f + trans[6]/Iscale + trans[7]/Iscale*(2.0f*norm/Iscale2-1.0f);
			float dfactor1 = 1.0f/Iscale;
			float dfactor2 = 1.0f/Iscale*(2.0f*norm/Iscale2-1.0f);
			
			dX[0][0] = ( factor*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( factor*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( factor*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = factor/rtx/scale;
			dX[0][6] = ( dfactor1*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			dX[0][7] = ( dfactor2*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][0] = ( factor*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( factor*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( factor*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = factor/rty/scale;
			dX[1][6] = ( dfactor1*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			dX[1][7] = ( dfactor2*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][0] = ( factor*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( factor*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( factor*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = factor/rtz/scale;
			dX[2][6] = ( dfactor1*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
			dX[2][7] = ( dfactor2*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==TRI_QUADRATIC_SCALE) {
			float normx = rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3]; float normx2 = normx*normx;
			float normy = rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4]; float normy2 = normy*normy;
			float normz = rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5]; float normz2 = normz*normz;

			float factorx = 1.0f + trans[6]/Iscale + trans[ 9]/Iscale*(2.0f*normx2/Iscale2-1.0f);
			float factory = 1.0f + trans[7]/Iscale + trans[10]/Iscale*(2.0f*normy2/Iscale2-1.0f);
			float factorz = 1.0f + trans[8]/Iscale + trans[11]/Iscale*(2.0f*normz2/Iscale2-1.0f);
						
			float dfactor1 = 1.0f/Iscale;
						
			float dfactorx2 = 1.0f/Iscale*(2.0f*normx2/Iscale2-1.0f);
			float dfactory2 = 1.0f/Iscale*(2.0f*normy2/Iscale2-1.0f);
			float dfactorz2 = 1.0f/Iscale*(2.0f*normz2/Iscale2-1.0f);
						
			dX[0][ 0] = ( ( factorx + 4.0f*trans[ 9]/Iscale*normx2/Iscale2 )*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][ 1] = ( ( factorx + 4.0f*trans[ 9]/Iscale*normx2/Iscale2 )*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][ 2] = ( ( factorx + 4.0f*trans[ 9]/Iscale*normx2/Iscale2 )*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][ 3] = ( ( factorx + 4.0f*trans[ 9]/Iscale*normx2/Iscale2 )/rtx )/scale;
			dX[0][ 6] = (  dfactor1*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			dX[0][ 9] = ( dfactorx2*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][ 0] = ( ( factory + 4.0f*trans[10]/Iscale*normy2/Iscale2 )*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][ 1] = ( ( factory + 4.0f*trans[10]/Iscale*normy2/Iscale2 )*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][ 2] = ( ( factory + 4.0f*trans[10]/Iscale*normy2/Iscale2 )*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][ 4] = ( ( factory + 4.0f*trans[10]/Iscale*normy2/Iscale2 )/rty )/scale;
			dX[1][ 7] = (  dfactor1*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			dX[1][10] = ( dfactory2*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][ 0] = ( ( factorz + 4.0f*trans[11]/Iscale*normz2/Iscale2 )*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][ 1] = ( ( factorz + 4.0f*trans[11]/Iscale*normz2/Iscale2 )*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][ 2] = ( ( factorz + 4.0f*trans[11]/Iscale*normz2/Iscale2 )*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][ 5] = ( ( factorz + 4.0f*trans[11]/Iscale*normz2/Iscale2 )/rtz )/scale;
			dX[2][ 8] = (  dfactor1*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
			dX[2][11] = ( dfactorz2*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==FULLY_AFFINE) {
			dX[0][0] = ( Xi/rtx )/scale;
			dX[0][1] = ( Yi/rtx )/scale;
			dX[0][2] = ( Zi/rtx )/scale;
			dX[0][9] = 1/rtx/scale;
			
			dX[1][3] = ( Xi/rty )/scale;
			dX[1][4] = ( Yi/rty )/scale;
			dX[1][5] = ( Zi/rty )/scale;
			dX[1][10] = 1/rty/scale;
			
			dX[2][6] = ( Xi/rtz )/scale;
			dX[2][7] = ( Yi/rtz )/scale;
			dX[2][8] = ( Zi/rtz )/scale;
			dX[2][11] = 1/rtz/scale;
		} else if (type==FULLY_QUADRATIC) {
			dX[0][0] = ( Xi/rtx )/scale;
			dX[0][1] = ( Yi/rtx )/scale;
			dX[0][2] = ( Zi/rtx )/scale;
			dX[0][9] = 1.0f/rtx/scale;
			dX[0][12] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rtx )/scale;
			dX[0][13] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rtx )/scale;
			dX[0][14] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rtx )/scale;
			dX[0][15] = ( Xi*Yi/Iscale/rtx )/scale;
			dX[0][16] = ( Yi*Zi/Iscale/rtx )/scale;
			dX[0][17] = ( Zi*Xi/Iscale/rtx )/scale;
			
			dX[1][3] = ( Xi/rty )/scale;
			dX[1][4] = ( Yi/rty )/scale;
			dX[1][5] = ( Zi/rty )/scale;
			dX[1][10] = 1.0f/rty/scale;
			dX[1][18] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rty )/scale;
			dX[1][19] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rty )/scale;
			dX[1][20] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rty )/scale;
			dX[1][21] = ( Xi*Yi/Iscale/rty )/scale;
			dX[1][22] = ( Yi*Zi/Iscale/rty )/scale;
			dX[1][23] = ( Zi*Xi/Iscale/rty )/scale;
			
			dX[2][6] = ( Xi/rtz )/scale;
			dX[2][7] = ( Yi/rtz )/scale;
			dX[2][8] = ( Zi/rtz )/scale;
			dX[2][11] = 1.0f/rtz/scale;
			dX[2][24] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rtz )/scale;
			dX[2][25] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rtz )/scale;
			dX[2][26] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rtz )/scale;
			dX[2][27] = ( Xi*Yi/Iscale/rtz )/scale;
			dX[2][28] = ( Yi*Zi/Iscale/rtz )/scale;
			dX[2][29] = ( Zi*Xi/Iscale/rtz )/scale;
		} else if (type==FULLY_CUBIC) {
			dX[0][0] = ( Xi/rtx )/scale;
			dX[0][1] = ( Yi/rtx )/scale;
			dX[0][2] = ( Zi/rtx )/scale;
			dX[0][9] = 1.0f/rtx/scale;
			dX[0][12] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rtx )/scale;
			dX[0][13] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rtx )/scale;
			dX[0][14] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rtx )/scale;
			dX[0][15] = ( Xi*Yi/Iscale/rtx )/scale;
			dX[0][16] = ( Yi*Zi/Iscale/rtx )/scale;
			dX[0][17] = ( Zi*Xi/Iscale/rtx )/scale;
			dX[0][30] = ( (4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi)/rtx )/scale;
			dX[0][31] = ( (4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi)/rtx )/scale;
			dX[0][32] = ( (4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)/rtx )/scale;
			dX[0][33] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale/rtx )/scale;
			dX[0][34] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale/rtx )/scale;
			dX[0][35] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale/rtx )/scale;
			dX[0][36] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale/rtx )/scale;
			dX[0][37] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale/rtx )/scale;
			dX[0][38] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale/rtx )/scale;
			dX[0][39] = ( Xi*Yi*Zi/Iscale2/rtx )/scale;
			
			dX[1][3] = ( Xi/rty )/scale;
			dX[1][4] = ( Yi/rty )/scale;
			dX[1][5] = ( Zi/rty )/scale;
			dX[1][10] = 1.0f/rty/scale;
			dX[1][18] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rty )/scale;
			dX[1][19] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rty )/scale;
			dX[1][20] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rty )/scale;
			dX[1][21] = ( Xi*Yi/Iscale/rty )/scale;
			dX[1][22] = ( Yi*Zi/Iscale/rty )/scale;
			dX[1][23] = ( Zi*Xi/Iscale/rty )/scale;
			dX[1][40] = ( (4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi)/rty )/scale;
			dX[1][41] = ( (4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi)/rty )/scale;
			dX[1][42] = ( (4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)/rty )/scale;
			dX[1][43] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale/rty )/scale;
			dX[1][44] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale/rty )/scale;
			dX[1][45] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale/rty )/scale;
			dX[1][46] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale/rty )/scale;
			dX[1][47] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale/rty )/scale;
			dX[1][48] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale/rty )/scale;
			dX[1][49] = ( Xi*Yi*Zi/Iscale2/rty )/scale;
			
			dX[2][6] = ( Xi/rtz )/scale;
			dX[2][7] = ( Yi/rtz )/scale;
			dX[2][8] = ( Zi/rtz )/scale;
			dX[2][11] = 1.0f/rtz/scale;
			dX[2][24] = ( (2.0f*Xi*Xi/Iscale-Iscale)/rtz )/scale;
			dX[2][25] = ( (2.0f*Yi*Yi/Iscale-Iscale)/rtz )/scale;
			dX[2][26] = ( (2.0f*Zi*Zi/Iscale-Iscale)/rtz )/scale;
			dX[2][27] = ( Xi*Yi/Iscale/rtz )/scale;
			dX[2][28] = ( Yi*Zi/Iscale/rtz )/scale;
			dX[2][29] = ( Zi*Xi/Iscale/rtz )/scale;
			dX[1][50] = ( (4.0f*Xi*Xi*Xi/Iscale2-3.0f*Xi)/rtz )/scale;
			dX[1][51] = ( (4.0f*Yi*Yi*Yi/Iscale2-3.0f*Yi)/rtz )/scale;
			dX[1][52] = ( (4.0f*Zi*Zi*Zi/Iscale2-3.0f*Zi)/rtz )/scale;
			dX[1][53] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Yi/Iscale/rtz )/scale;
			dX[1][54] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Zi/Iscale/rtz )/scale;
			dX[1][55] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Xi/Iscale/rtz )/scale;
			dX[1][56] = ( (2.0f*Xi*Xi/Iscale-Iscale)*Zi/Iscale/rtz )/scale;
			dX[1][57] = ( (2.0f*Yi*Yi/Iscale-Iscale)*Xi/Iscale/rtz )/scale;
			dX[1][58] = ( (2.0f*Zi*Zi/Iscale-Iscale)*Yi/Iscale/rtz )/scale;
			dX[1][59] = ( Xi*Yi*Zi/Iscale2/rtz )/scale;
			
		}
		return dX;
	}

	/** 
	 *	computes the transformed coordinates from image to template space
	 *	with scaling (2: half res, 4: quarter res, etc)
	 */
	public final float normScalingFactor(int x,int y,int z, float[] trans, float[][] rot, float scale) {
		float scaling = 1.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
			scaling = 1.0f;
		} else if (type==RIGID) {
			scaling = 1.0f;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			scaling = factor*factor*factor;
		} else if (type==SCALED_RIGID) {
			scaling = (1.0f+trans[6]/Iscale)*(1.0f+trans[7]/Iscale)*(1.0f+trans[8]/Iscale);
		} else if (type==QUADRATIC_SCALE) {
			float norm = Xi*Xi + Yi*Yi + Zi*Zi;
			float factor = 1.0f + trans[6]/Iscale + trans[7]/Iscale*(2.0f*norm/Iscale2-1.0f);
			scaling = factor*factor*factor;
		} else if (type==TRI_QUADRATIC_SCALE) {
			float normx = rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3]; normx *= normx;
			float normy = rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4]; normy *= normy;
			float normz = rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5]; normz *= normz;

			float factorx = 1.0f + trans[6]/Iscale + trans[ 9]/Iscale*(2.0f*normx/Iscale2-1.0f);
			float factory = 1.0f + trans[7]/Iscale + trans[10]/Iscale*(2.0f*normy/Iscale2-1.0f);
			float factorz = 1.0f + trans[8]/Iscale + trans[11]/Iscale*(2.0f*normz/Iscale2-1.0f);
			scaling = factorx*factory*factorz;
		} else if (type==FULLY_AFFINE) {
			scaling = (1+trans[0])*(1+trans[4])*(1+trans[8])
						+ trans[1]*trans[5]*trans[6]
						+ trans[2]*trans[3]*trans[7]
						- (1+trans[0])*trans[5]*trans[7]
						- trans[2]*(1+trans[4])*trans[6]
						- trans[1]*trans[3]*(1+trans[8]);
		} else if (type==FULLY_QUADRATIC) {
			scaling = 1.0f;
		} else if (type==FULLY_CUBIC) {
			scaling = 1.0f;
		}
		return scaling;
	}

	public final float[] normScalingFactorDerivatives(int x,int y,int z,float[] trans, float[][] rot, float[][] dRa, float[][] dRb, float[][] dRc, float scale) {
		float[] dF = new float[Nt];
		
		for (int t=0;t<Nt;t++)
			dF[t] = 0.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
		} else if (type==RIGID) {
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dF[6] = 3.0f*dfactor*factor*factor;
		} else if (type==SCALED_RIGID) {
			float factorx = 1.0f + trans[6]/Iscale;
			float factory = 1.0f + trans[7]/Iscale;
			float factorz = 1.0f + trans[8]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dF[6] = dfactor*factory*factorz;
			dF[7] = dfactor*factorz*factorx;
			dF[8] = dfactor*factorx*factory;
		} else if (type==QUADRATIC_SCALE) {
			float norm = Xi*Xi + Yi*Yi + Zi*Zi;
			float factor = 1.0f + trans[6]/Iscale + trans[7]/Iscale*(2.0f*norm/Iscale2-1.0f);
			float dfactor1 = 1.0f/Iscale;
			float dfactor2 = 1.0f/Iscale*(2.0f*norm/Iscale2-1.0f);
			
			dF[6] = 3.0f*dfactor1*factor*factor;
			dF[7] = 3.0f*dfactor2*factor*factor;
		} else if (type==TRI_QUADRATIC_SCALE) {
			float normx = rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3]; float normx2 = normx*normx;
			float normy = rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4]; float normy2 = normy*normy;
			float normz = rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5]; float normz2 = normz*normz;

			float factorx = 1.0f + trans[6]/Iscale + trans[ 9]/Iscale*(2.0f*normx2/Iscale2-1.0f);
			float factory = 1.0f + trans[7]/Iscale + trans[10]/Iscale*(2.0f*normy2/Iscale2-1.0f);
			float factorz = 1.0f + trans[8]/Iscale + trans[11]/Iscale*(2.0f*normz2/Iscale2-1.0f);
						
			float dfactor1 = 1.0f/Iscale;
						
			float dfactorx2 = 1.0f/Iscale*(2.0f*normx2/Iscale2-1.0f);
			float dfactory2 = 1.0f/Iscale*(2.0f*normy2/Iscale2-1.0f);
			float dfactorz2 = 1.0f/Iscale*(2.0f*normz2/Iscale2-1.0f);
						
			dF[ 6] =  dfactor1*factory*factorz;
			dF[ 9] = dfactorx2*factory*factorz;
			
			dF[ 7] =  dfactor1*factorz*factorx;
			dF[10] = dfactory2*factorz*factorx;
			
			dF[ 8] =  dfactor1*factorx*factory;
			dF[11] = dfactorz2*factorx*factory;
		} else if (type==FULLY_AFFINE) {
			dF[0] = (1+trans[4])*(1+trans[8]) - trans[5]*trans[7];
			dF[1] = trans[5]*trans[6] - trans[3]*(1+trans[8]);
			dF[2] = trans[3]*trans[7] - (1+trans[4])*trans[6];
			dF[3] = trans[2]*trans[7] - trans[1]*(1+trans[8]);
			dF[4] = (1+trans[0])*(1+trans[8]) - trans[2]*trans[6];
			dF[5] = trans[1]*trans[6] - (1+trans[0])*trans[7];
			dF[6] = trans[1]*trans[5] - trans[2]*(1+trans[4]);
			dF[7] = trans[2]*trans[3] - (1+trans[0])*trans[5];
			dF[8] = (1+trans[0])*(1+trans[4]) - trans[1]*trans[3];
		} else if (type==FULLY_QUADRATIC) {
		} else if (type==FULLY_CUBIC) {
		}
		return dF;
	}

	public final void precomputeImageToTemplateMatrix(float[][] transformMatrix, float[] trans, float[][] rot, float scale) {
		
		if (type==NONE) {
			transformMatrix[0][0] = rix/rtx;
			transformMatrix[0][1] = 0;
			transformMatrix[0][2] = 0;
			transformMatrix[1][0] = 0;
			transformMatrix[1][1] = riy/rty;
			transformMatrix[1][2] = 0;
			transformMatrix[2][0] = 0;
			transformMatrix[2][1] = 0;
			transformMatrix[2][2] = riz/rtz;
			
			transformMatrix[0][3] = ((-x0i)*rix/rtx + x0t)/scale;
			transformMatrix[1][3] = ((-y0i)*riy/rty + y0t)/scale;
			transformMatrix[2][3] = ((-z0i)*riz/rtz + z0t)/scale;		
		} else if (type==RIGID) {
			transformMatrix[0][0] = rot[0][0]*rix/rtx;
			transformMatrix[0][1] = rot[0][1]*riy/rtx;
			transformMatrix[0][2] = rot[0][2]*riz/rtx;
			transformMatrix[1][0] = rot[1][0]*rix/rty;
			transformMatrix[1][1] = rot[1][1]*riy/rty;
			transformMatrix[1][2] = rot[1][2]*riz/rty;
			transformMatrix[2][0] = rot[2][0]*rix/rtz;
			transformMatrix[2][1] = rot[2][1]*riy/rtz;
			transformMatrix[2][2] = rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = ((rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = ((rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = ((rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f+trans[6]/Iscale;
			
			transformMatrix[0][0] = factor*rot[0][0]*rix/rtx;
			transformMatrix[0][1] = factor*rot[0][1]*riy/rtx;
			transformMatrix[0][2] = factor*rot[0][2]*riz/rtx;
			transformMatrix[1][0] = factor*rot[1][0]*rix/rty;
			transformMatrix[1][1] = factor*rot[1][1]*riy/rty;
			transformMatrix[1][2] = factor*rot[1][2]*riz/rty;
			transformMatrix[2][0] = factor*rot[2][0]*rix/rtz;
			transformMatrix[2][1] = factor*rot[2][1]*riy/rtz;
			transformMatrix[2][2] = factor*rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = (factor*(rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = (factor*(rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = (factor*(rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		}  else if (type==SCALED_RIGID) {
			float factorx = 1.0f+trans[6]/Iscale;
			float factory = 1.0f+trans[7]/Iscale;
			float factorz = 1.0f+trans[8]/Iscale;
			
			transformMatrix[0][0] = factorx*rot[0][0]*rix/rtx;
			transformMatrix[0][1] = factorx*rot[0][1]*riy/rtx;
			transformMatrix[0][2] = factorx*rot[0][2]*riz/rtx;
			transformMatrix[1][0] = factory*rot[1][0]*rix/rty;
			transformMatrix[1][1] = factory*rot[1][1]*riy/rty;
			transformMatrix[1][2] = factory*rot[1][2]*riz/rty;
			transformMatrix[2][0] = factorz*rot[2][0]*rix/rtz;
			transformMatrix[2][1] = factorz*rot[2][1]*riy/rtz;
			transformMatrix[2][2] = factorz*rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = (factorx*(rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = (factory*(rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = (factorz*(rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		} else if (type==FULLY_AFFINE) {
			transformMatrix[0][0] = (1+trans[0])*rix/rtx;
			transformMatrix[0][1] = (0+trans[1])*riy/rtx;
			transformMatrix[0][2] = (0+trans[2])*riz/rtx;
			transformMatrix[1][0] = (0+trans[3])*rix/rty;
			transformMatrix[1][1] = (1+trans[4])*riy/rty;
			transformMatrix[1][2] = (0+trans[5])*riz/rty;
			transformMatrix[2][0] = (0+trans[6])*rix/rtz;
			transformMatrix[2][1] = (0+trans[7])*riy/rtz;
			transformMatrix[2][2] = (1+trans[8])*riz/rtz;
			
			transformMatrix[0][3] = ((((1+trans[0])*(-x0i)*rix+   trans[1] *(-y0i)*riy+   trans[2] *(-z0i)*riz)+trans[ 9])/rtx + x0t)/scale;
			transformMatrix[1][3] = (((   trans[3] *(-x0i)*rix+(1+trans[4])*(-y0i)*riy+   trans[5] *(-z0i)*riz)+trans[10])/rty + y0t)/scale;
			transformMatrix[2][3] = (((   trans[6] *(-x0i)*rix+   trans[7] *(-y0i)*riy+(1+trans[8])*(-z0i)*riz)+trans[11])/rtz + z0t)/scale;
		} 	
		
		// debug info
		/*
		MedicUtil.displayMessage("transform to matrix :\n - original transform: "+displayTransform(trans)
								+"img center ("+x0i+", "+y0i+", "+z0i+")\n"
								+"tpl center ("+x0t+", "+y0t+", "+z0t+")\n");
		*/
		return;
	}
	
	public final RotationMatrix computeRotationMatrix(float[] trans) {
		RotationMatrix R = new RotationMatrix();
		R.setParameters(trans[0],trans[1],trans[2]);
		return R;
	}
	
	public static final float[][] computeRotation(float[] trans) {
		RotationMatrix R = new RotationMatrix();
		R.setParameters(trans[0],trans[1],trans[2]);
		return R.getMatrix();
	}
	
	public static final float[] changeTransformType(float[] trans, String oldType_, String newType_) {
		int oldType, newType;
		int No, Nn;
		float[] newtrans=null;
		
			 if (oldType_.equals("rigid")) { oldType = RIGID; No = 6; } 
		else if (oldType_.equals("single_scale")) { oldType = SINGLE_SCALE; No = 7; }
		else if (oldType_.equals("scaled_rigid")) { oldType = SCALED_RIGID; No = 9; }
		else if (oldType_.equals("fully_affine")) { oldType = FULLY_AFFINE; No = 12; }
		else if (oldType_.equals("quadratic_scale")) { oldType = QUADRATIC_SCALE; No = 8; }
		else if (oldType_.equals("tri_quadratic_scale")) { oldType = TRI_QUADRATIC_SCALE; No = 12; }
		else if (oldType_.equals("fully_quadratic")) { oldType = FULLY_QUADRATIC; No = 30; }
		else if (oldType_.equals("fully_cubic")) { oldType = FULLY_CUBIC; No = 60; }
		else { oldType = NONE; No = 0; }

			 if (newType_.equals("rigid")) { newType = RIGID; Nn = 6; } 
		else if (newType_.equals("single_scale")) { newType = SINGLE_SCALE; Nn = 7; }
		else if (newType_.equals("scaled_rigid")) { newType = SCALED_RIGID; Nn = 9; }
		else if (newType_.equals("fully_affine")) { newType = FULLY_AFFINE; Nn = 12; }
		else if (newType_.equals("quadratic_scale")) { newType = QUADRATIC_SCALE; Nn = 8; }
		else if (newType_.equals("tri_quadratic_scale")) { newType = TRI_QUADRATIC_SCALE; Nn = 12; }
		else if (newType_.equals("fully_quadratic")) { newType = FULLY_QUADRATIC; Nn = 30; }
		else if (newType_.equals("fully_cubic")) { newType = FULLY_CUBIC; Nn = 60; }
		else { newType = NONE; Nn = 0; }

		if (oldType==newType) return trans;
		else if (oldType==NONE) {
			newtrans = new float[Nn];
			for (int n=0;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==RIGID && (newType==SINGLE_SCALE || newType==SCALED_RIGID || newType==QUADRATIC_SCALE || newType==TRI_QUADRATIC_SCALE) ) {
			newtrans = new float[Nn];
			for (int n=0;n<No;n++) newtrans[n] = trans[n];
			for (int n=No;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==SINGLE_SCALE && newType==QUADRATIC_SCALE) {
			newtrans = new float[Nn];
			for (int n=0;n<No;n++) newtrans[n] = trans[n];
			for (int n=No;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==SINGLE_SCALE && (newType==SCALED_RIGID || newType==TRI_QUADRATIC_SCALE) ) {
			newtrans = new float[Nn];
			for (int n=0;n<No;n++) newtrans[n] = trans[n];
			newtrans[No] = trans[No-1];
			newtrans[No+1] = trans[No-1];
			for (int n=No+2;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==RIGID && (newType==FULLY_AFFINE || newType==FULLY_QUADRATIC || newType==FULLY_CUBIC) ) {
			float[][] rot = computeRotation(trans);
			newtrans = new float[Nn];
			for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
				if (n==m) newtrans[n+3*m] = rot[m][n]-1.0f;
				else newtrans[n+3*m] = rot[m][n];
			}
			for (int n=0;n<3;n++) newtrans[9+n] = trans[3+n];
			for (int n=12;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==FULLY_AFFINE && (newType==FULLY_QUADRATIC || newType==FULLY_CUBIC) ) {
			newtrans = new float[Nn];
			for (int n=0;n<12;n++) newtrans[n] = trans[n];
			for (int n=12;n<Nn;n++) newtrans[n] = 0.0f;
		} else {
			MedicUtilPublic.displayMessage("change from "+oldType_+" to "+newType_+" currently not supported\n");
			return trans;
		}
		MedicUtilPublic.displayMessage("change from "+oldType_+" to "+newType_+" done ("+newtrans.length+")\n");
		return newtrans;
	}
	
	/** 
	 *	computes transformed directions from template to image space
	 */
	public final void templateToImageDirection(float[] dir, int x,int y,int z, float[] trans, float[][] rot, float scale) {
		float norm = (float)Math.sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
		
		if (norm<1e-30f) return;
		
		float[] V = new float[3];
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		if (type==NONE) {
			return;
		} else if (type==RIGID) {
			V[0] = (rot[0][0]*dir[0]+rot[1][0]*dir[1]+rot[2][0]*dir[2]);
			V[1] = (rot[0][1]*dir[0]+rot[1][1]*dir[1]+rot[2][1]*dir[2]);
			V[2] = (rot[0][2]*dir[0]+rot[1][2]*dir[1]+rot[2][2]*dir[2]);
		} else if (type==SINGLE_SCALE) {
			V[0] = (rot[0][0]*dir[0]+rot[1][0]*dir[1]+rot[2][0]*dir[2]);
			V[1] = (rot[0][1]*dir[0]+rot[1][1]*dir[1]+rot[2][1]*dir[2]);
			V[2] = (rot[0][2]*dir[0]+rot[1][2]*dir[1]+rot[2][2]*dir[2]);
		} else if (type==SCALED_RIGID) {
			float sx = 1.0f/(1.0f+trans[6]/Iscale);
			float sy = 1.0f/(1.0f+trans[7]/Iscale);
			float sz = 1.0f/(1.0f+trans[8]/Iscale);
			
			V[0] = (rot[0][0]*dir[0]/sx+rot[1][0]*dir[1]/sy+rot[2][0]*dir[2]/sz);
			V[1] = (rot[0][1]*dir[0]/sx+rot[1][1]*dir[1]/sy+rot[2][1]*dir[2]/sz);
			V[2] = (rot[0][2]*dir[0]/sx+rot[1][2]*dir[1]/sy+rot[2][2]*dir[2]/sz);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		} else if (type==QUADRATIC_SCALE) {
			V[0] = (rot[0][0]*dir[0]+rot[1][0]*dir[1]+rot[2][0]*dir[2]);
			V[1] = (rot[0][1]*dir[0]+rot[1][1]*dir[1]+rot[2][1]*dir[2]);
			V[2] = (rot[0][2]*dir[0]+rot[1][2]*dir[1]+rot[2][2]*dir[2]);
		} else if (type==TRI_QUADRATIC_SCALE) {
			float normx = rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3]; normx *= normx;
			float normy = rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4]; normy *= normy;
			float normz = rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5]; normz *= normz;

			float sx = 1.0f + trans[6]/Iscale + trans[ 9]/Iscale*(2.0f*normx/Iscale2-1.0f);
			float sy = 1.0f + trans[7]/Iscale + trans[10]/Iscale*(2.0f*normy/Iscale2-1.0f);
			float sz = 1.0f + trans[8]/Iscale + trans[11]/Iscale*(2.0f*normz/Iscale2-1.0f);
			
			V[0] = (rot[0][0]*dir[0]/sx+rot[1][0]*dir[1]/sy+rot[2][0]*dir[2]/sz);
			V[1] = (rot[0][1]*dir[0]/sx+rot[1][1]*dir[1]/sy+rot[2][1]*dir[2]/sz);
			V[2] = (rot[0][2]*dir[0]/sx+rot[1][2]*dir[1]/sy+rot[2][2]*dir[2]/sz);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		} else if (type==FULLY_AFFINE) {
			float[][] M = new float[3][3];
			M[0][0] = (1+trans[0]);
			M[0][1] = trans[1];
			M[0][2] = trans[2];
			M[1][0] = trans[3];
			M[1][1] = (1+trans[4]);
			M[1][2] = trans[5];
			M[2][0] = trans[6];
			M[2][1] = trans[7];
			M[2][2] = (1+trans[8]);
			Numerics.invert3Dmatrix(M);
			
			V[0] = (M[0][0]*dir[0]+M[0][1]*dir[1]+M[0][2]*dir[2]);
			V[1] = (M[1][0]*dir[0]+M[1][1]*dir[1]+M[1][2]*dir[2]);
			V[2] = (M[2][0]*dir[0]+M[2][1]*dir[1]+M[2][2]*dir[2]);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
			
		} else if (type==FULLY_QUADRATIC) {
			float[][] M = new float[3][3];
			M[0][0] = (1+trans[0]) +trans[12]*2.0f*Xi/Iscale +0.5f*trans[15]*Yi/Iscale +0.5f*trans[17]*Zi/Iscale;
			M[0][1] =    trans[1]  +trans[13]*2.0f*Yi/Iscale +0.5f*trans[15]*Xi/Iscale +0.5f*trans[16]*Zi/Iscale;
			M[0][2] =    trans[2]  +trans[14]*2.0f*Zi/Iscale +0.5f*trans[16]*Yi/Iscale +0.5f*trans[17]*Xi/Iscale;
			M[1][0] =    trans[3]  +trans[18]*2.0f*Xi/Iscale +0.5f*trans[21]*Yi/Iscale +0.5f*trans[23]*Zi/Iscale;
			M[1][1] = (1+trans[4]) +trans[19]*2.0f*Yi/Iscale +0.5f*trans[21]*Xi/Iscale +0.5f*trans[22]*Zi/Iscale;
			M[1][2] =    trans[5]  +trans[20]*2.0f*Zi/Iscale +0.5f*trans[22]*Yi/Iscale +0.5f*trans[23]*Xi/Iscale;
			M[2][0] =    trans[6]  +trans[24]*2.0f*Xi/Iscale +0.5f*trans[27]*Yi/Iscale +0.5f*trans[29]*Zi/Iscale;
			M[2][1] =    trans[7]  +trans[25]*2.0f*Yi/Iscale +0.5f*trans[27]*Xi/Iscale +0.5f*trans[28]*Zi/Iscale;
			M[2][2] = (1+trans[8]) +trans[26]*2.0f*Zi/Iscale +0.5f*trans[28]*Yi/Iscale +0.5f*trans[29]*Xi/Iscale;
			Numerics.invert3Dmatrix(M);
			
			V[0] = (M[0][0]*dir[0]+M[0][1]*dir[1]+M[0][2]*dir[2]);
			V[1] = (M[1][0]*dir[0]+M[1][1]*dir[1]+M[1][2]*dir[2]);
			V[2] = (M[2][0]*dir[0]+M[2][1]*dir[1]+M[2][2]*dir[2]);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		} else if (type==FULLY_CUBIC) {
			float[][] M = new float[3][3];
			M[0][0] = (1+trans[0]) +trans[12]*2.0f*Xi/Iscale +0.5f*trans[15]*Yi/Iscale +0.5f*trans[17]*Zi/Iscale +trans[30]*(4.0f*Xi*Xi/Iscale2-3.0f) +0.6667f*trans[33]*2.0f*Xi*Yi/Iscale2 +0.3333f*trans[35]*2.0f*Zi*Zi/Iscale2 +0.6667f*trans[36]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[37]*2.0f*Yi*Yi/Iscale2 +0.3333f*trans[39]*Yi*Zi/Iscale2;
			M[0][1] =    trans[1]  +trans[13]*2.0f*Yi/Iscale +0.5f*trans[15]*Xi/Iscale +0.5f*trans[16]*Zi/Iscale +trans[31]*(4.0f*Yi*Yi/Iscale2-3.0f) +0.3333f*trans[33]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[34]*2.0f*Yi*Zi/Iscale2 +0.6667f*trans[37]*2.0f*Yi*Xi/Iscale2 +0.3333f*trans[38]*2.0f*Zi*Zi/Iscale2 +0.3333f*trans[39]*Xi*Zi/Iscale2;
			M[0][2] =    trans[2]  +trans[14]*2.0f*Zi/Iscale +0.5f*trans[16]*Yi/Iscale +0.5f*trans[17]*Xi/Iscale +trans[32]*(4.0f*Zi*Zi/Iscale2-3.0f) +0.3333f*trans[34]*2.0f*Yi*Yi/Iscale2 +0.6667f*trans[35]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[36]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[38]*2.0f*Yi*Zi/Iscale2 +0.3333f*trans[39]*Xi*Yi/Iscale2;
			M[1][0] =    trans[3]  +trans[18]*2.0f*Xi/Iscale +0.5f*trans[21]*Yi/Iscale +0.5f*trans[23]*Zi/Iscale +trans[40]*(4.0f*Xi*Xi/Iscale2-3.0f) +0.6667f*trans[43]*2.0f*Xi*Yi/Iscale2 +0.3333f*trans[45]*2.0f*Zi*Zi/Iscale2 +0.6667f*trans[46]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[47]*2.0f*Yi*Yi/Iscale2 +0.3333f*trans[49]*Yi*Zi/Iscale2;
			M[1][1] = (1+trans[4]) +trans[19]*2.0f*Yi/Iscale +0.5f*trans[21]*Xi/Iscale +0.5f*trans[22]*Zi/Iscale +trans[41]*(4.0f*Yi*Yi/Iscale2-3.0f) +0.3333f*trans[43]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[44]*2.0f*Yi*Zi/Iscale2 +0.6667f*trans[47]*2.0f*Yi*Xi/Iscale2 +0.3333f*trans[48]*2.0f*Zi*Zi/Iscale2 +0.3333f*trans[49]*Xi*Zi/Iscale2;
			M[1][2] =    trans[5]  +trans[20]*2.0f*Zi/Iscale +0.5f*trans[22]*Yi/Iscale +0.5f*trans[23]*Xi/Iscale +trans[42]*(4.0f*Zi*Zi/Iscale2-3.0f) +0.3333f*trans[44]*2.0f*Yi*Yi/Iscale2 +0.6667f*trans[45]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[46]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[48]*2.0f*Yi*Zi/Iscale2 +0.3333f*trans[49]*Xi*Yi/Iscale2;
			M[2][0] =    trans[6]  +trans[24]*2.0f*Xi/Iscale +0.5f*trans[27]*Yi/Iscale +0.5f*trans[29]*Zi/Iscale +trans[50]*(4.0f*Xi*Xi/Iscale2-3.0f) +0.6667f*trans[53]*2.0f*Xi*Yi/Iscale2 +0.3333f*trans[55]*2.0f*Zi*Zi/Iscale2 +0.6667f*trans[56]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[57]*2.0f*Yi*Yi/Iscale2 +0.3333f*trans[59]*Yi*Zi/Iscale2;
			M[2][1] =    trans[7]  +trans[25]*2.0f*Yi/Iscale +0.5f*trans[27]*Xi/Iscale +0.5f*trans[28]*Zi/Iscale +trans[51]*(4.0f*Yi*Yi/Iscale2-3.0f) +0.3333f*trans[53]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[54]*2.0f*Yi*Zi/Iscale2 +0.6667f*trans[57]*2.0f*Yi*Xi/Iscale2 +0.3333f*trans[58]*2.0f*Zi*Zi/Iscale2 +0.3333f*trans[59]*Xi*Zi/Iscale2;
			M[2][2] = (1+trans[8]) +trans[26]*2.0f*Zi/Iscale +0.5f*trans[28]*Yi/Iscale +0.5f*trans[29]*Xi/Iscale +trans[52]*(4.0f*Zi*Zi/Iscale2-3.0f) +0.3333f*trans[54]*2.0f*Yi*Yi/Iscale2 +0.6667f*trans[55]*2.0f*Xi*Zi/Iscale2 +0.3333f*trans[56]*2.0f*Xi*Xi/Iscale2 +0.6667f*trans[58]*2.0f*Yi*Zi/Iscale2 +0.3333f*trans[59]*Xi*Yi/Iscale2;
			Numerics.invert3Dmatrix(M);
			
			V[0] = (M[0][0]*dir[0]+M[0][1]*dir[1]+M[0][2]*dir[2]);
			V[1] = (M[1][0]*dir[0]+M[1][1]*dir[1]+M[1][2]*dir[2]);
			V[2] = (M[2][0]*dir[0]+M[2][1]*dir[1]+M[2][2]*dir[2]);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		} 
		dir[0] = V[0];
		dir[1] = V[1];
		dir[2] = V[2];
		
		return;
	}

	public final String displayTransform(float[] trans) {
		String info = "transform: (";
		for (int n=0;n<Nt-1;n++) info += trans[n]+", ";
		info += trans[Nt-1]+")\n";
		
		return info;
	}
	
}
