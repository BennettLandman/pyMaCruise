package edu.jhu.ece.iacl.algorithms.gvf;

import java.util.List;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector3f;

import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;
import edu.jhu.ece.iacl.algorithms.tgdm.GenericTGDM.NBPoint;
import edu.jhu.ece.iacl.algorithms.tgdm.TrackingTGDM.EmbeddedNBPoint;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;

/**
 * Compute regular gradient or upwind gradient for a volumetric function.
 * 
 * @author Blake Lucas
 *
 */
public class FastMarchingGradient extends AbstractCalculation {
	public enum Normalization {
		NONE, MAGNITUDE, MAG_SQUARED
	};

	public FastMarchingGradient() {
		super();
		setLabel("Fast-Marching Gradient");
	}

	public FastMarchingGradient(AbstractCalculation parent) {
		super(parent);
		setLabel("Fast-Marching Gradient");
	}

	// Compute normal velocity between two level sets
	public float[][][][] velocity(float[][][] levelMat1, float[][][] levelMat2,
			float dt) {
		float[][][][] grad = gradient(levelMat1, Normalization.MAG_SQUARED);
		int rows = levelMat1.length;
		int cols = levelMat1[0].length;
		int slices = levelMat1[0][0].length;
		float dpdt;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					dpdt = (levelMat1[i][j][k] - levelMat2[i][j][k]) / dt;
					grad[i][j][k][0] *= dpdt;
					grad[i][j][k][1] *= dpdt;
					grad[i][j][k][2] *= dpdt;
				}
			}
		}
		return grad;
	}
	public static Vector3f interpolate(Point3f p,float[][][][] gradient){
		int i0, j0, k0, i1, j1, k1;
		double x=p.x;
		double y=p.y;
		double z=p.z;
		Vector3f v=new Vector3f();
		int rows=gradient.length;
		int cols=gradient[0].length;
		int slices=gradient[0][0].length;
		double dx, dy, dz, hx, hy, hz;
		j0 = (int) Math.floor(x);
		i0 = (int) Math.floor(y);
		k0 = (int) Math.floor(z);
		j1 = j0+1;
		i1 = i0+1;
		k1 = k0+1;
		if (j0 < 0 || j1 > (rows - 1) || i0 < 0 || i1 > (cols - 1) || k0 < 0
				|| k1 > (slices - 1)) {
			System.out.println("INVALID "+x+" "+y+" "+z+" "+rows+" "+cols+" "+slices);
			j0=(int)Math.min(Math.max(j0, 0),rows-1);
			i0=(int)Math.min(Math.max(i0, 0),cols-1);
			k0=(int)Math.min(Math.max(k0, 0),slices-1);

			float[] vec=new float[3];
			v.x=gradient[j0][i0][k0][0];
			v.y=gradient[j0][i0][k0][1];
			v.z=gradient[j0][i0][k0][2];
			return v;
		} else {
			dx = x - j0;
			dy = y - i0;
			dz = z - k0;
			/* Introduce more variables to reduce computation */
			hx = 1.0 - dx;
			hy = 1.0 - dy;
			hz = 1.0 - dz;
			v.x=(float)(((gradient[j0][i0][k0][0] * hx + 
					gradient[j1][i0][k0][0]* dx)* hy + 
					(gradient[j0][i1][k0][0]* hx + 
							gradient[j1][i1][k0][0]* dx)* dy)* hz+ 
					((gradient[j0][i0][k1][0] * hx + 
							gradient[j1][i0][k1][0]* dx)* hy + 
					(gradient[j0][i1][k1][0] * hx + 
							gradient[j1][i1][k1][0]* dx)
							* dy) * dz);
			v.y=(float)(((gradient[j0][i0][k0][1] * hx + 
					gradient[j1][i0][k0][1]* dx)* hy + 
					(gradient[j0][i1][k0][1]* hx + 
							gradient[j1][i1][k0][1]* dx)* dy)* hz+ 
					((gradient[j0][i0][k1][1] * hx + 
							gradient[j1][i0][k1][1]* dx)* hy + 
					(gradient[j0][i1][k1][1] * hx + 
							gradient[j1][i1][k1][1]* dx)
							* dy) * dz);
			v.z=(float)(((gradient[j0][i0][k0][2] * hx + 
					gradient[j1][i0][k0][2]* dx)* hy + 
					(gradient[j0][i1][k0][2]* hx + 
							gradient[j1][i1][k0][2]* dx)* dy)* hz+ 
					((gradient[j0][i0][k1][2] * hx + 
							gradient[j1][i0][k1][2]* dx)* hy + 
					(gradient[j0][i1][k1][2] * hx + 
							gradient[j1][i1][k1][2]* dx)
							* dy) * dz);
			
			double mag=Math.sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
			if(mag>0){
				v.x/=mag;
				v.y/=mag;
				v.z/=mag;
			} else {
				System.out.println("MAGNITUDE ZERO "+x+","+y+","+z);
			}
			
			return v;
		}
}	
	public float[][][][] gradient(float[][][] levelMat, Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Dx, Dy, Dz;
		float[][][][] gradMat = new float[rows][cols][slices][3];
		Point3f p = new Point3f();
		setTotalUnits(rows);
		for (int i = 0; i < rows; i++) {
			incrementCompletedUnits();
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					LX = (i == 0) ? 1 : 0;
					HX = (i == (rows - 1)) ? 1 : 0;

					LY = (j == 0) ? 1 : 0;
					HY = (j == (cols - 1)) ? 1 : 0;

					LZ = (k == 0) ? 1 : 0;
					HZ = (k == (slices - 1)) ? 1 : 0;

					Nv = levelMat[i][j - 1 + LY][k];
					Sv = levelMat[i][j + 1 - HY][k];
					Wv = levelMat[i - 1 + LX][j][k];
					Ev = levelMat[i + 1 - HX][j][k];
					Fv = levelMat[i][j][k + 1 - HZ];
					Bv = levelMat[i][j][k - 1 + LZ];
					Cv = levelMat[i][j][k];

					if (Cv > 0) {

						p.x = (Wv < Cv) ? ((Wv < Ev) ? Cv - Wv : Ev - Cv)
								: ((Ev < Cv) ? Ev - Cv : 0);
						p.y = (Nv < Cv) ? ((Nv < Sv) ? Cv - Nv : Sv - Cv)
								: ((Sv < Cv) ? Sv - Cv : 0);
						p.z = (Bv < Cv) ? ((Bv < Fv) ? Cv - Bv : Fv - Cv)
								: ((Fv < Cv) ? Fv - Cv : 0);
					} else {
						p.x = (Wv > Cv) ? ((Wv > Ev) ? Cv - Wv : Ev - Cv)
								: ((Ev > Cv) ? Ev - Cv : 0);
						p.y = (Nv > Cv) ? ((Nv > Sv) ? Cv - Nv : Sv - Cv)
								: ((Sv > Cv) ? Sv - Cv : 0);
						p.z = (Bv > Cv) ? ((Bv > Fv) ? Cv - Bv : Fv - Cv)
								: ((Fv > Cv) ? Fv - Cv : 0);
					}
					switch (norm) {
					case NONE:
						break;
					case MAGNITUDE:
						GeometricUtilities.normalize(p);
						break;
					case MAG_SQUARED:
						GeometricUtilities.normalizeSqr(p);
						break;
					}
					gradMat[i][j][k][0] = p.x;
					gradMat[i][j][k][1] = p.y;
					gradMat[i][j][k][2] = p.z;
				}
			}
		}
		markCompleted();
		return gradMat;
	}

	public void velocity(float[][][] lastLevel, float[][][] nextLevel,
			List<EmbeddedNBPoint> NBANDList) {
		int rows = lastLevel.length;
		int cols = lastLevel[0].length;
		int slices = lastLevel[0][0].length;
		float dpdt;
		int i, j, k;
		int n = 0;
		FastMarchingGradient fmg=new FastMarchingGradient(this);
		float[][] grad=fmg.gradient(lastLevel, NBANDList,Normalization.MAG_SQUARED);
		for (EmbeddedNBPoint pt : NBANDList) {
			i = pt.x;
			j = pt.y;
			k = pt.z;
			dpdt = (lastLevel[i][j][k] - nextLevel[i][j][k]);
			pt.vec=new Vector3d(grad[n][0],grad[n][1],grad[n][2]);
			pt.vec.x *= dpdt;
			pt.vec.y *= dpdt;
			pt.vec.z *= dpdt;
			n++;
		}
	}
	public Vector3d gradientUpwind(float[][][] levelMat,
			float[][][] targetMat, NBPoint pt,
			Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Nq, Sq, Wq, Eq, Fq, Bq, Cq;
		float Dx, Dy, Dz;
		Vector3d grad = new Vector3d();
		int i, j, k;
		int n = 0;
		double maxVal = 0;

			i = pt.x;
			j = pt.y;
			k = pt.z;
			LX = (i == 0) ? 1 : 0;
			HX = (i == (rows - 1)) ? 1 : 0;

			LY = (j == 0) ? 1 : 0;
			HY = (j == (cols - 1)) ? 1 : 0;

			LZ = (k == 0) ? 1 : 0;
			HZ = (k == (slices - 1)) ? 1 : 0;

			Nv = levelMat[i][j - 1 + LY][k];
			Sv = levelMat[i][j + 1 - HY][k];
			Wv = levelMat[i - 1 + LX][j][k];
			Ev = levelMat[i + 1 - HX][j][k];
			Fv = levelMat[i][j][k + 1 - HZ];
			Bv = levelMat[i][j][k - 1 + LZ];
			Cv = levelMat[i][j][k];

			Nq = targetMat[i][j - 1 + LY][k];
			Sq = targetMat[i][j + 1 - HY][k];
			Wq = targetMat[i - 1 + LX][j][k];
			Eq = targetMat[i + 1 - HX][j][k];
			Fq = targetMat[i][j][k + 1 - HZ];
			Bq = targetMat[i][j][k - 1 + LZ];
			Cq = targetMat[i][j][k];
			if(Float.isNaN(Cq)){
				grad.x=0;
				grad.y=0;
				grad.z=0;
			} else {
				if (Float.isNaN(Nq))
					Nv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Sq))
					Sv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Eq))
					Ev = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Wq))
					Wv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Fq))
					Fv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Bq))
					Bv = Math.signum(Cv) * Float.MAX_VALUE;
	
				/*
				 * if(Cv>0){ p.x=(((Wq>=-1&&Wv<Cv)?Cq-Wq:0)+((Eq>=-1&&Ev<Cv)?Eq-Cq:0))/(1+((Wq>=-1&&Eq>=-1&&Wv<Cv&&Ev<Cv)?1:0));
				 * p.y=(((Nq>=-1&&Nv<Cv)?Cq-Nq:0)+((Sq>=-1&&Sv<Cv)?Sq-Cq:0))/(1+((Nq>=-1&&Sq>=-1&&Nv<Cv&&Sv<Cv)?1:0));
				 * p.z=(((Bq>=-1&&Bv<Cv)?Cq-Bq:0)+((Fq>=-1&&Fv<Cv)?Fq-Cq:0))/(1+((Bq>=-1&&Fq>=-1&&Fv<Cv&&Bv<Cv)?1:0)); }
				 * else {
				 * p.x=(((Wq>=-1&&Wv>Cv)?Cq-Wq:0)+((Eq>=-1&&Ev>Cv)?Eq-Cq:0))/(1+((Wq>=-1&&Eq>=-1&&Wv>Cv&&Ev>Cv)?1:0));
				 * p.y=(((Nq>=-1&&Nv>Cv)?Cq-Nq:0)+((Sq>=-1&&Sv>Cv)?Sq-Cq:0))/(1+((Nq>=-1&&Sq>=-1&&Nv>Cv&&Sv>Cv)?1:0));
				 * p.z=(((Bq>=-1&&Bv>Cv)?Cq-Bq:0)+((Fq>=-1&&Fv>Cv)?Fq-Cq:0))/(1+((Bq>=-1&&Fq>=-1&&Fv>Cv&&Bv>Cv)?1:0)); }
				 */
				if (Cv > 0) {
	
					grad.x = (Wv < Cv) ? ((Wv < Ev) ? Cq - Wq : Eq - Cq)
							: ((Ev < Cv) ? Eq - Cq : 0);
					grad.y = (Nv < Cv) ? ((Nv < Sv) ? Cq - Nq : Sq - Cq)
							: ((Sv < Cv) ? Sq - Cq : 0);
					grad.z = (Bv < Cv) ? ((Bv < Fv) ? Cq - Bq : Fq - Cq)
							: ((Fv < Cv) ? Fq - Cq : 0);
				} else {
					grad.x = (Wv > Cv) ? ((Wv > Ev) ? Cq - Wq : Cq - Eq)
							: ((Ev > Cv) ? Cq - Eq : 0);
					grad.y = (Nv > Cv) ? ((Nv > Sv) ? Cq - Nq : Sq - Cq)
							: ((Sv > Cv) ? Sq - Cq : 0);
					grad.z = (Bv > Cv) ? ((Bv > Fv) ? Cq - Bq : Fq - Cq)
							: ((Fv > Cv) ? Fq - Cq : 0);
				}
			}
			/*
			 * if(GeomUtil.length(p)>2&&Cv>0){
			 * System.out.println(GeomUtil.length(p)+" Cv "+Cv+" Nv "+Nv+" Sv
			 * "+Sv+" Wv "+Wv+" Ev "+Ev+" Fv "+Fv+" Bv "+Bv);
			 * System.out.println("Cq "+Cq+" Nq "+Nq+" Sq "+Sq+" Wq "+Wq+" Eq
			 * "+Eq+" Fq "+Fq+" Bq "+Bq); }
			 */
			switch (norm) {
			case NONE:
				break;
			case MAGNITUDE:
				GeometricUtilities.normalize(grad);
				break;
			case MAG_SQUARED:
				GeometricUtilities.normalizeSqr(grad);
				break;
			}
		return grad;
	}
	public float[][] gradientUpwind(float[][][] levelMat,
			float[][][] targetMat, List<NBPoint> NBANDList,
			Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Nq, Sq, Wq, Eq, Fq, Bq, Cq;
		float Dx, Dy, Dz;
		float[][] gradMat = new float[NBANDList.size()][3];
		Point3f p = new Point3f();
		int i, j, k;
		int n = 0;
		double maxVal = 0;
		for (NBPoint pt : NBANDList) {
			i = pt.x;
			j = pt.y;
			k = pt.z;
			LX = (i == 0) ? 1 : 0;
			HX = (i == (rows - 1)) ? 1 : 0;

			LY = (j == 0) ? 1 : 0;
			HY = (j == (cols - 1)) ? 1 : 0;

			LZ = (k == 0) ? 1 : 0;
			HZ = (k == (slices - 1)) ? 1 : 0;

			Nv = levelMat[i][j - 1 + LY][k];
			Sv = levelMat[i][j + 1 - HY][k];
			Wv = levelMat[i - 1 + LX][j][k];
			Ev = levelMat[i + 1 - HX][j][k];
			Fv = levelMat[i][j][k + 1 - HZ];
			Bv = levelMat[i][j][k - 1 + LZ];
			Cv = levelMat[i][j][k];

			Nq = targetMat[i][j - 1 + LY][k];
			Sq = targetMat[i][j + 1 - HY][k];
			Wq = targetMat[i - 1 + LX][j][k];
			Eq = targetMat[i + 1 - HX][j][k];
			Fq = targetMat[i][j][k + 1 - HZ];
			Bq = targetMat[i][j][k - 1 + LZ];
			Cq = targetMat[i][j][k];
			if(Float.isNaN(Cq)){
				p.x=0;
				p.y=0;
				p.z=0;
			} else {
				if (Float.isNaN(Nq))
					Nv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Sq))
					Sv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Eq))
					Ev = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Wq))
					Wv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Fq))
					Fv = Math.signum(Cv) * Float.MAX_VALUE;
				if (Float.isNaN(Bq))
					Bv = Math.signum(Cv) * Float.MAX_VALUE;
	
				/*
				 * if(Cv>0){ p.x=(((Wq>=-1&&Wv<Cv)?Cq-Wq:0)+((Eq>=-1&&Ev<Cv)?Eq-Cq:0))/(1+((Wq>=-1&&Eq>=-1&&Wv<Cv&&Ev<Cv)?1:0));
				 * p.y=(((Nq>=-1&&Nv<Cv)?Cq-Nq:0)+((Sq>=-1&&Sv<Cv)?Sq-Cq:0))/(1+((Nq>=-1&&Sq>=-1&&Nv<Cv&&Sv<Cv)?1:0));
				 * p.z=(((Bq>=-1&&Bv<Cv)?Cq-Bq:0)+((Fq>=-1&&Fv<Cv)?Fq-Cq:0))/(1+((Bq>=-1&&Fq>=-1&&Fv<Cv&&Bv<Cv)?1:0)); }
				 * else {
				 * p.x=(((Wq>=-1&&Wv>Cv)?Cq-Wq:0)+((Eq>=-1&&Ev>Cv)?Eq-Cq:0))/(1+((Wq>=-1&&Eq>=-1&&Wv>Cv&&Ev>Cv)?1:0));
				 * p.y=(((Nq>=-1&&Nv>Cv)?Cq-Nq:0)+((Sq>=-1&&Sv>Cv)?Sq-Cq:0))/(1+((Nq>=-1&&Sq>=-1&&Nv>Cv&&Sv>Cv)?1:0));
				 * p.z=(((Bq>=-1&&Bv>Cv)?Cq-Bq:0)+((Fq>=-1&&Fv>Cv)?Fq-Cq:0))/(1+((Bq>=-1&&Fq>=-1&&Fv>Cv&&Bv>Cv)?1:0)); }
				 */
				if (Cv > 0) {
	
					p.x = (Wv < Cv) ? ((Wv < Ev) ? Cq - Wq : Eq - Cq)
							: ((Ev < Cv) ? Eq - Cq : 0);
					p.y = (Nv < Cv) ? ((Nv < Sv) ? Cq - Nq : Sq - Cq)
							: ((Sv < Cv) ? Sq - Cq : 0);
					p.z = (Bv < Cv) ? ((Bv < Fv) ? Cq - Bq : Fq - Cq)
							: ((Fv < Cv) ? Fq - Cq : 0);
				} else {
					p.x = (Wv > Cv) ? ((Wv > Ev) ? Cq - Wq : Cq - Eq)
							: ((Ev > Cv) ? Cq - Eq : 0);
					p.y = (Nv > Cv) ? ((Nv > Sv) ? Cq - Nq : Sq - Cq)
							: ((Sv > Cv) ? Sq - Cq : 0);
					p.z = (Bv > Cv) ? ((Bv > Fv) ? Cq - Bq : Fq - Cq)
							: ((Fv > Cv) ? Fq - Cq : 0);
				}
			}
			/*
			 * if(GeomUtil.length(p)>2&&Cv>0){
			 * System.out.println(GeomUtil.length(p)+" Cv "+Cv+" Nv "+Nv+" Sv
			 * "+Sv+" Wv "+Wv+" Ev "+Ev+" Fv "+Fv+" Bv "+Bv);
			 * System.out.println("Cq "+Cq+" Nq "+Nq+" Sq "+Sq+" Wq "+Wq+" Eq
			 * "+Eq+" Fq "+Fq+" Bq "+Bq); }
			 */
			switch (norm) {
			case NONE:
				break;
			case MAGNITUDE:
				GeometricUtilities.normalize(p);
				break;
			case MAG_SQUARED:
				GeometricUtilities.normalizeSqr(p);
				break;
			}
			gradMat[n][0] = p.x;
			gradMat[n][1] = p.y;
			gradMat[n][2] = p.z;
			n++;
		}
		// System.out.println("MAX VAL "+maxVal);
		return gradMat;
	}
	public float[][] gradientCentral(float[][][] levelMat, List<NBPoint> NBANDList,
			Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Nq, Sq, Wq, Eq, Fq, Bq, Cq;
		float Dx, Dy, Dz;
		float[][] gradMat = new float[NBANDList.size()][3];
		Point3f p = new Point3f();
		int i, j, k;
		int n = 0;
		double maxVal = 0;
		for (NBPoint pt : NBANDList) {
			i = pt.x;
			j = pt.y;
			k = pt.z;
			LX = (i == 0) ? 1 : 0;
			HX = (i == (rows - 1)) ? 1 : 0;

			LY = (j == 0) ? 1 : 0;
			HY = (j == (cols - 1)) ? 1 : 0;

			LZ = (k == 0) ? 1 : 0;
			HZ = (k == (slices - 1)) ? 1 : 0;

			Nv = levelMat[i][j - 1 + LY][k];
			Sv = levelMat[i][j + 1 - HY][k];
			Wv = levelMat[i - 1 + LX][j][k];
			Ev = levelMat[i + 1 - HX][j][k];
			Fv = levelMat[i][j][k + 1 - HZ];
			Bv = levelMat[i][j][k - 1 + LZ];
			Cv = levelMat[i][j][k];

			if(Float.isNaN(Cv)||Float.isNaN(Nv)||Float.isNaN(Sv)||Float.isNaN(Wv)||Float.isNaN(Ev)||Float.isNaN(Fv)||Float.isNaN(Bv)){
				p.x=0;
				p.y=0;
				p.z=0;
			} else {
				p.x = 0.5f*(Ev-Wv);
				p.y = 0.5f*(Sv-Nv);
				p.z = 0.5f*(Fv-Bv);
			}
			switch (norm) {
			case NONE:
				break;
			case MAGNITUDE:
				GeometricUtilities.normalize(p);
				break;
			case MAG_SQUARED:
				GeometricUtilities.normalizeSqr(p);
				break;
			}
			gradMat[n][0] = p.x;
			gradMat[n][1] = p.y;
			gradMat[n][2] = p.z;
			n++;
		}
		// System.out.println("MAX VAL "+maxVal);
		return gradMat;
	}

	public float[][][][] gradientUpwind(float[][][] levelMat,
			float[][][] targetMat, Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Nq, Sq, Wq, Eq, Fq, Bq, Cq;
		float Dx, Dy, Dz;
		float[][][][] gradMat = new float[rows][cols][slices][3];
		Point3f p = new Point3f();
		setTotalUnits(rows);
		for (int i = 0; i < rows; i++) {
			incrementCompletedUnits();
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {

					LX = (i == 0) ? 1 : 0;
					HX = (i == (rows - 1)) ? 1 : 0;

					LY = (j == 0) ? 1 : 0;
					HY = (j == (cols - 1)) ? 1 : 0;

					LZ = (k == 0) ? 1 : 0;
					HZ = (k == (slices - 1)) ? 1 : 0;

					Nv = levelMat[i][j - 1 + LY][k];
					Sv = levelMat[i][j + 1 - HY][k];
					Wv = levelMat[i - 1 + LX][j][k];
					Ev = levelMat[i + 1 - HX][j][k];
					Fv = levelMat[i][j][k + 1 - HZ];
					Bv = levelMat[i][j][k - 1 + LZ];
					Cv = levelMat[i][j][k];

					Nq = targetMat[i][j - 1 + LY][k];
					Sq = targetMat[i][j + 1 - HY][k];
					Wq = targetMat[i - 1 + LX][j][k];
					Eq = targetMat[i + 1 - HX][j][k];
					Fq = targetMat[i][j][k + 1 - HZ];
					Bq = targetMat[i][j][k - 1 + LZ];
					Cq = targetMat[i][j][k];
					if (Cv > 0) {

						p.x = (Wv < Cv) ? ((Wv < Ev) ? Cq - Wq : Eq - Cq)
								: ((Ev < Cv) ? Eq - Cq : 0);
						p.y = (Nv < Cv) ? ((Nv < Sv) ? Cq - Nq : Sq - Cq)
								: ((Sv < Cv) ? Sq - Cq : 0);
						p.z = (Bv < Cv) ? ((Bv < Fv) ? Cq - Bq : Fq - Cq)
								: ((Fv < Cv) ? Fq - Cq : 0);
					} else {
						p.x = (Wv > Cv) ? ((Wv > Ev) ? Cq - Wq : Eq - Cq)
								: ((Ev > Cv) ? Eq - Cq : 0);
						p.y = (Nv > Cv) ? ((Nv > Sv) ? Cq - Nq : Sq - Cq)
								: ((Sv > Cv) ? Sq - Cq : 0);
						p.z = (Bv > Cv) ? ((Bv > Fv) ? Cq - Bq : Fq - Cq)
								: ((Fv > Cv) ? Fq - Cq : 0);
					}
					switch (norm) {
					case NONE:
						break;
					case MAGNITUDE:
						GeometricUtilities.normalize(p);
						break;
					case MAG_SQUARED:
						GeometricUtilities.normalizeSqr(p);
						break;
					}
					gradMat[i][j][k][0] = p.x;
					gradMat[i][j][k][1] = p.y;
					gradMat[i][j][k][2] = p.z;
				}
			}
		}
		return gradMat;
	}

	public float[][] gradient(float[][][] levelMat,
			List<EmbeddedNBPoint> NBANDList, Normalization norm) {
		int rows = levelMat.length;
		int cols = levelMat[0].length;
		int slices = levelMat[0][0].length;
		int LX, HX, LY, HY, LZ, HZ;
		float Nv, Sv, Wv, Ev, Fv, Bv, Cv;
		float Dx, Dy, Dz;
		float[][] gradMat = new float[NBANDList.size()][3];
		Point3f p = new Point3f();
		int i, j, k;
		int n = 0;
		double minG=1E30,maxG=-1E30;
		for (NBPoint pt : NBANDList) {
			i = pt.x;
			j = pt.y;
			k = pt.z;
			LX = (i == 0) ? 1 : 0;
			HX = (i == (rows - 1)) ? 1 : 0;

			LY = (j == 0) ? 1 : 0;
			HY = (j == (cols - 1)) ? 1 : 0;

			LZ = (k == 0) ? 1 : 0;
			HZ = (k == (slices - 1)) ? 1 : 0;

			Nv = levelMat[i][j - 1 + LY][k];
			Sv = levelMat[i][j + 1 - HY][k];
			Wv = levelMat[i - 1 + LX][j][k];
			Ev = levelMat[i + 1 - HX][j][k];
			Fv = levelMat[i][j][k + 1 - HZ];
			Bv = levelMat[i][j][k - 1 + LZ];
			Cv = levelMat[i][j][k];

			if (Cv > 0) {

				p.x = (Wv < Cv) ? ((Wv < Ev) ? Cv - Wv : Ev - Cv)
						: ((Ev < Cv) ? Ev - Cv : 0);
				p.y = (Nv < Cv) ? ((Nv < Sv) ? Cv - Nv : Sv - Cv)
						: ((Sv < Cv) ? Sv - Cv : 0);
				p.z = (Bv < Cv) ? ((Bv < Fv) ? Cv - Bv : Fv - Cv)
						: ((Fv < Cv) ? Fv - Cv : 0);
			} else {
				p.x = (Wv > Cv) ? ((Wv > Ev) ? Cv - Wv : Ev - Cv)
						: ((Ev > Cv) ? Ev - Cv : 0);
				p.y = (Nv > Cv) ? ((Nv > Sv) ? Cv - Nv : Sv - Cv)
						: ((Sv > Cv) ? Sv - Cv : 0);
				p.z = (Bv > Cv) ? ((Bv > Fv) ? Cv - Bv : Fv - Cv)
						: ((Fv > Cv) ? Fv - Cv : 0);
			}
			double l=GeometricUtilities.length(p);
			if(l!=0){
				minG=Math.min(minG,l);
				maxG=Math.max(maxG,l);
			}
			switch (norm) {
				case NONE:
					break;
				case MAGNITUDE:
					GeometricUtilities.normalize(p);
					break;
				case MAG_SQUARED:
					GeometricUtilities.normalizeSqr(p);
					break;
			}
			gradMat[n][0] = p.x;
			gradMat[n][1] = p.y;
			gradMat[n][2] = p.z;
			
			n++;
		}
		System.out.println("MIN GRAD "+minG+" MAX GRAD "+maxG);
		return gradMat;
	}
}
