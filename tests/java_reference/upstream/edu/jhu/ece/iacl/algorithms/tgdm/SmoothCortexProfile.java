package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;

public class SmoothCortexProfile extends ProfileTGDM{
	protected double levelSet;
	public SmoothCortexProfile(double curvatureForce,double externalForce,double pressureForce,double levelSet,int maxIters){
		
		super(4.7f,3.6f,4.7f+1.5f,5.0f,30,30,0.1f,0.01f,curvatureForce,externalForce,pressureForce,maxIters);
		this.levelSet=levelSet;
		this.skipLastIters=false;
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
		 return curvatureForce*K*GPhi;
	}
	@Override
	public Vector3d getEffectiveVelocity(double K, double G,
			Vector3d curvatureGrad, double regionForce, Vector3d regionGrad,
			Vector3d advectionForce, Vector3d advectionGrad) {
		double GPhi=curvatureGrad.length();
		double tmp;
		Vector3d vel;
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
			  vel=new Vector3d(curvatureGrad);
			  GeometricUtilities.normalize(vel);
			  vel.scale(curvatureForce*K);
		} else {
			K = 0; /* 0/0 */
			vel=new Vector3d();
		}
		return vel;
	}
}
