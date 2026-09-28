package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;

public class DynamicPressureCortexProfile extends ProfileTGDM{

	public DynamicPressureCortexProfile(double curvatureForce,double externalForce,double pressureForce,int maxIters){
		super(4.8f,3.6f,4.8f+1.2f,1000f,30,30,0.1f,0.045f,curvatureForce,externalForce,pressureForce,maxIters);	
		finalMaxSignDist=5;
		// TODO Auto-generated constructor stub
	}
	public double getEffectiveForce(double K,double G,Vector3d curvatureGrad,double regionForce,Vector3d regionGrad,Vector3d advectionForce,Vector3d advectionGrad){
		double GPhi=curvatureGrad.length();
		GPhi = GPhi + 0.001f;	
		double tmp = GPhi * GPhi;
		tmp = GPhi*GPhi;
		K = K/(tmp*GPhi); 
		/* Note: this is 2 times the true mean curvature*/
		G = G/(tmp*tmp);
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		 return curvatureForce*K*GPhi+regionForce*regionGrad.length();
	}
	public Vector3d getEffectiveVelocity(double K, double G,
			Vector3d curvatureGrad, double regionForce, Vector3d regionGrad,
			Vector3d advectionForce, Vector3d advectionGrad) {
		Vector3d vel;
		double GPhi=curvatureGrad.length();
		if(GPhi > 0.00001){
			vel=new Vector3d(curvatureGrad);
			vel.normalize();
		} else {
			vel=new Vector3d();
		}
		GPhi = GPhi + 0.001f;	
		double tmp = GPhi * GPhi;
		tmp = GPhi*GPhi;
		K = K/(tmp*GPhi); 
		/* Note: this is 2 times the true mean curvature*/
		G = G/(tmp*tmp);
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		vel.scale(curvatureForce*K);
		Vector3d regionV=new Vector3d(regionGrad);
		GeometricUtilities.normalize(regionV);
		regionV.scale(regionForce);
		vel.add(regionV);
		vel.add(advectionForce);
		return vel;
	}
/*
	public double getLevelSetChange(double phi,double cu,double cv,double cw,double cg, double GPhi, double K, double G,double Vp,double Dmx,double Dpx,double Dmy,double Dpy,double Dmz,double Dpz) {
		GPhi = GPhi + 0.001f;	
		double tmp = GPhi * GPhi;
		double value,DeltaP,DeltaM;
		tmp = GPhi*GPhi;
		K = K/(tmp*GPhi); 
		// Note: this is 2 times the true mean curvature
		G = G/(tmp*tmp);
		
		// Now K is the weighted curvatue 
		tmp = K*K - 2*G;
		if(tmp > 0)
		  K = K*G/tmp;
		
		// Compute gradient in downwind direction 
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

		// Compute gradient in upwind direction 
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
		value = (Vp > 0)?DeltaP:DeltaM;
		value = value*Vp*cg + curvatureForce*K*GPhi; // Vp is already Vp*Kp 
		// value = cg*value; 
		
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
