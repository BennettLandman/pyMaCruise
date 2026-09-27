package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;

public class CentralCortexProfile extends ProfileTGDM{
	public CentralCortexProfile(double curvatureForce,double externalForce,double pressureForce,int maxIters){
		super(4.7f,3.6f,4.7f+1.5f,12f,30,10,0.1f,0.02f,curvatureForce,externalForce,pressureForce,maxIters);	
		// TODO Auto-generated constructor stub
	}
	public double getEffectiveForce(double K,double G,Vector3d curvatureGrad,double regionForce,Vector3d regionGrad,Vector3d advectionForce,Vector3d advectionGrad){
		double GPhi=curvatureGrad.length();
		double tmp;
		if(GPhi > 0.00001){ 
			  tmp = GPhi*GPhi;
			  K = K/(tmp*GPhi); 
			  /* Note: this is 2 times the true mean curvature*/
			  G = G/(tmp*tmp);
			  /* printf("K = %g, G = %g \n", K, G); */
			  /* Now K is the weighted curvatue */
			  tmp = K*K - 2*G;
			  if(tmp > 0)
			    K = K*G/tmp;
		} else K = 0; /* 0/0 */
		 return curvatureForce*K*GPhi+regionForce*regionGrad.length()+advectionForce.dot(advectionGrad);
	}
	public Vector3d getEffectiveVelocity(double K, double G,
			Vector3d curvatureGrad, double regionForce, Vector3d regionGrad,
			Vector3d advectionForce, Vector3d advectionGrad) {
		Vector3d vel;
		double GPhi=curvatureGrad.length();
		double tmp;
		//Augment phi gradient
		if(GPhi > 0.00001){ 
			  vel=new Vector3d(curvatureGrad);
			  vel.normalize();
			  tmp = GPhi*GPhi;
			  K = K/(tmp*GPhi); 
			  /* Note: this is 2 times the true mean curvature*/
			  G = G/(tmp*tmp);
			  /* printf("K = %g, G = %g \n", K, G); */
			  /* Now K is the weighted curvatue */
			  tmp = K*K - 2*G;
			  if(tmp > 0)
			    K = K*G/tmp;
		} else {
			K = 0; /* 0/0 */
			vel=new Vector3d();
		}
		GeometricUtilities.normalize(vel);
		vel.scale(curvatureForce*K);
		Vector3d regionV=new Vector3d(regionGrad);
		GeometricUtilities.normalize(regionV);
		regionV.scale(regionForce);
		vel.add(regionV);
		vel.add(advectionForce);
		return vel;
	}
	/*
	public double getLevelSetChange(double phi,double cu, double cv, double cw, double cg, double GPhi, double K, double G, double Vp, double Dmx, double Dpx, double Dmy, double Dpy, double Dmz, double Dpz) {
		double tmp,DeltaP,DeltaM,value;
		if(GPhi > 0.00001){ 
			  tmp = GPhi*GPhi;
			  K = K/(tmp*GPhi); 
			  G = G/(tmp*tmp);

			  tmp = K*K - 2*G;
			  if(tmp > 0)
			    K = K*G/tmp;
		} else K = 0; 
		// Compute DeltaP 
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
		
		// Add Pressure and Curvature terms 
		if(Vp > 0) value = DeltaP;
		else value = DeltaM;
		value = cg*(value*Vp + curvatureForce*K*GPhi); // Vp is already Vp*Kp 
		
		// Add External force term 
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
		
		return value;
	}
*/
}
