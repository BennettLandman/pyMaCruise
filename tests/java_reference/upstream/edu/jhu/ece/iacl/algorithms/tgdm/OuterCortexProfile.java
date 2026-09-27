package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;

public class OuterCortexProfile extends ProfileTGDM{

	public OuterCortexProfile(double curvatureForce,double externalForce,double pressureForce,int maxIters){
		super(4.7f,3.6f,4.7f+1.5f,15f,30,30,0.1f,0.045f,curvatureForce,externalForce,pressureForce,maxIters);	
		// TODO Auto-generated constructor stub
	}
	public double getEffectiveForce(double K,double G,Vector3d curvatureGrad,double regionForce,Vector3d regionGrad,Vector3d advectionForce,Vector3d advectionGrad){
		double GPhi=curvatureGrad.length();
		double tmp;
		double Ks;
		//Augment phi gradient
		GPhi = GPhi + 0.001f;	
		if(K < -1.5 && G < 0){
			regionForce = 0;
	        Ks = 3*curvatureForce;
		}else{
		  Ks = curvatureForce;
		}
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		if(K > 1.5) K = 1.5f;
		else if(K < -1.5) K = -1.5f;

		 return Ks*K*GPhi+regionForce*regionGrad.length()+advectionForce.dot(advectionGrad);
	}
	public Vector3d getEffectiveVelocity(double K, double G,
			Vector3d curvatureGrad, double regionForce, Vector3d regionGrad,
			Vector3d advectionForce, Vector3d advectionGrad) {
		Vector3d vel;
		double GPhi=curvatureGrad.length();
		double tmp;
		double Ks;
		//Augment phi gradient
		if(GPhi > 0.00001){
			vel=new Vector3d(curvatureGrad);
			vel.normalize();
		} else {
			vel=new Vector3d();
		}
		GPhi = GPhi + 0.001f;	
		if(K < -1.5 && G < 0){
			regionForce = 0;
	        Ks = 3*curvatureForce;
		}else{
		  Ks = curvatureForce;
		}
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		if(K > 1.5) K = 1.5f;
		else if(K < -1.5) K = -1.5f;
		vel.scale(Ks*K);
		Vector3d regionV=new Vector3d(regionGrad);
		GeometricUtilities.normalize(regionV);
		regionV.scale(regionForce);
		vel.add(regionV);
		vel.add(advectionForce);
		return vel;
	}
	/*
	public double getLevelSetChange(double phi,double cu, double cv, double cw, double cg, double GPhi, double K, double G, double Vp, double Dmx, double Dpx, double Dmy, double Dpy, double Dmz, double Dpz) {
		double tmp,DeltaP,DeltaM, value;
		double Ks;
		//Augment phi gradient
		GPhi = GPhi + 0.001f;	
		if(K < -1.5 && G < 0){
			Vp = 0;
	        Ks = 3*curvatureForce;
		}else{
		  Ks = curvatureForce;
		}
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		if(K > 1.5) K = 1.5f;
		else if(K < -1.5) K = -1.5f;
		
		
		tmp = (Dmx > 0) ? Dmx : 0;
		DeltaP = tmp*tmp;
		tmp = (Dpx < 0) ? Dpx : 0;
		DeltaP = DeltaP + tmp*tmp;
		tmp = (Dmy > 0) ? Dmy : 0;
		DeltaP = DeltaP + tmp*tmp;
		tmp = (Dpy < 0) ? Dpy : 0;
		DeltaP = DeltaP + tmp*tmp;
		tmp = (Dmz > 0) ? Dmz : 0;
		DeltaP = DeltaP + tmp*tmp;
		tmp = (Dpz < 0) ? Dpz : 0;
		DeltaP = DeltaP + tmp*tmp;
		DeltaP = (double)Math.sqrt(DeltaP);

		// Compute DeltaM 
		tmp = (Dmx < 0) ? Dmx : 0;
		DeltaM = tmp*tmp;
		tmp = (Dpx > 0) ? Dpx : 0;
		DeltaM = DeltaM + tmp*tmp;
		tmp = (Dmy < 0) ? Dmy : 0;
		DeltaM = DeltaM + tmp*tmp;
		tmp = (Dpy > 0) ? Dpy : 0;
		DeltaM = DeltaM + tmp*tmp;
		tmp = (Dmz < 0) ? Dmz : 0;
		DeltaM = DeltaM + tmp*tmp;
		tmp = (Dpz > 0) ? Dpz : 0;
		DeltaM = DeltaM + tmp*tmp;
		DeltaM = (double)Math.sqrt(DeltaM);


		//Add Pressure and Curvature terms 
		if(Vp > 0) value = DeltaP;
		else value = DeltaM;
		value = value*Vp + Ks*K*GPhi; // Vp is already Vp*Kp 
		value = cg*value;
		
		// Add External force term 
		if(externalForce != 0){
		  cu = externalForce*cu;
		  cv = externalForce*cv;
		  cw = externalForce*cw;
		  if(cu > 0)
		    value = value + cu*Dmx;
		  else
		    value = value + cu*Dpx;
		  if(cv > 0)
		    value = value + cv*Dmy;
		  else
		    value = value + cv*Dpy;
		  if(cw > 0)
		    value = value + cw*Dmz;
		  else
		    value = value + cw*Dpz;
		}
		return value;
	}
*/
}
