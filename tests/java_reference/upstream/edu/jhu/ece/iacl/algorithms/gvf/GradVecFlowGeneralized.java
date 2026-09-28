package edu.jhu.ece.iacl.algorithms.gvf;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.structures.image.ImageDataMath;
/**
 * Generalized solution to Gradient Vector Flow
 * @author Blake Lucas
 *
 */
public class GradVecFlowGeneralized extends GradVecFlow{
	private static final GradVecFlowGeneralized gvf=new GradVecFlowGeneralized();
	public static ImageData doSolve(ImageData f,float lambda,boolean normalize,int iter){
		return gvf.solve(f,lambda,normalize,iter);
	}
	public GradVecFlowGeneralized(AbstractCalculation parent) {
		super(parent);
	}
	public GradVecFlowGeneralized() {
		super();

	}

	private ImageDataFloat process(ImageDataFloat Gf,ImageDataFloat SqrMagf,float lambda,int ITERS){
		int _z, _y, _x;
		float sqrmagf;
		float lambda2;
		int xn=Gf.getRows();
		int yn=Gf.getCols();
		int zn=Gf.getSlices();
		float[][][] SqrMagfmat=SqrMagf.toArray3d();
		float[][][][] Gfmat=Gf.toArray4d();
		ImageDataFloat V=new ImageDataFloat(xn,yn,zn,3);
		float Vmat[][][][]=V.toArray4d();
		ImageDataFloat U=new ImageDataFloat(xn,yn,zn);
		float Umat[][][]=U.toArray3d();
		ImageDataFloat f=new ImageDataFloat(xn,yn,zn);
		float fmat[][][]=f.toArray3d();
		ImageDataFloat h=new ImageDataFloat(xn,yn,zn);
		float hmat[][][]=h.toArray3d();
		ImageDataFloat g=new ImageDataFloat(xn,yn,zn);
		float gmat[][][]=g.toArray3d();
		lambda2 = lambda*lambda;

		
		
		for (_x=0; _x<xn; _x++){
			for (_y=0; _y<yn; _y++){
				for (_z=0; _z<zn; _z++){
					sqrmagf = SqrMagfmat[_x][_y][_z];
					gmat[_x][_y][_z] = (float)Math.exp(-sqrmagf/lambda2);
					hmat[_x][_y][_z] = 1 - gmat[_x][_y][_z];
					gmat[_x][_y][_z] /= 6.0;
					fmat[_x][_y][_z] = 0 - hmat[_x][_y][_z]*Gfmat[_x][_y][_z][0];
					Umat[_x][_y][_z] = Gfmat[_x][_y][_z][0];
				}
			}
		}
		incrementCompletedUnits();
		Poisson3d poisson3dX=new Poisson3d(this);
		poisson3dX.setLabel("Poisson3d X Dimension");
		U=poisson3dX.solve(U, f, g, h,ITERS);
		Umat=U.toArray3d();
		for (_x=0; _x<xn; _x++){
			for (_y=0; _y<yn; _y++){
				for (_z=0; _z<zn; _z++){
					Vmat[_x][_y][_z][0]=Umat[_x][_y][_z];
					fmat[_x][_y][_z] = 0 - hmat[_x][_y][_z]*Gfmat[_x][_y][_z][1];
					Umat[_x][_y][_z] = Gfmat[_x][_y][_z][1];
				}
			}
		}
		incrementCompletedUnits();
		Poisson3d poisson3dY=new Poisson3d(this);
		poisson3dY.setLabel("Poisson3d Y Dimension");
		U=poisson3dY.solve(U, f, g, h,ITERS);
		Umat=U.toArray3d();
		for (_x=0; _x<xn; _x++){
			for (_y=0; _y<yn; _y++){
				for (_z=0; _z<zn; _z++){
					Vmat[_x][_y][_z][1]=Umat[_x][_y][_z];
					fmat[_x][_y][_z] = 0- hmat[_x][_y][_z]*Gfmat[_x][_y][_z][2];
					Umat[_x][_y][_z] = Gfmat[_x][_y][_z][2];
				}
			}
		}
		incrementCompletedUnits();

		Poisson3d poisson3dZ=new Poisson3d(this);
		poisson3dZ.setLabel("Poisson3d Z Dimension");
		U=poisson3dZ.solve(U, f, g, h,ITERS);		
		Umat=U.toArray3d();

		for (_x=0; _x<xn; _x++){
			for (_y=0; _y<yn; _y++){
				for (_z=0; _z<zn; _z++){
					Vmat[_x][_y][_z][2]=Umat[_x][_y][_z];
				}
			}
		}
		markCompleted();
		return V;
	}
	/**
	 * Solve GVF using generalized method
	 * @param vol matrix 3d 
	 * @param lambda lambda
	 * @param normalize normalize if true
	 * @param iter maximum iterations
	 * @return vector matrix3d
	 */
	public ImageData solve(ImageData _f,float lambda,boolean normalize,int iter){
			ImageDataFloat f=new ImageDataFloat(_f);

			int xn=f.getRows();
			int yn=f.getCols();
			int zn=f.getSlices();
			System.out.println("Lambda "+lambda+" normalize "+normalize+" iters "+iter);
			
			//if(_f.getType()==VoxelType.UBYTE){
				float volmin=1000,volmax=0;
				float[][][] tmp=f.toArray3d();
				for(int i=0;i<xn;i++){
					for(int j=0;j<yn;j++){
						for(int k=0;k<zn;k++){
							volmin=Math.min(volmin,tmp[i][j][k]);
							volmax=Math.max(volmax,tmp[i][j][k]);
						}
					}
				}				
				System.out.println("Normalize Volume "+volmin+" "+volmax);
				
				for(int i=0;i<xn;i++){
					for(int j=0;j<yn;j++){
						for(int k=0;k<zn;k++){
							tmp[i][j][k]=(tmp[i][j][k]-volmin)/(volmax-volmin);
						}
					}
				}
			//}
			setLabel("Gradient Vector Flow");
			String name=f.getName()+"_gvf";
			ImageDataFloat V;//=new VectorMatrixX(xn,yn,zn,new Vector3());  
			ImageDataFloat Gf;//=new VectorMatrixX(xn,yn,zn,new Vector3());
			ImageDataFloat SqrMagf;//=new Matrix3(xn,yn,zn);
			setTotalUnits(5);
			Gf=Gradient3d.doSolve(f);   /* compute the gradient of ivol     */

			//calculated euclidian norm
//			SqrMagf=Gf.mag();
			
			SqrMagf=ImageDataMath.vectorMagnitude4D(Gf);
			incrementCompletedUnits();
			V=process(Gf,SqrMagf,lambda,iter);
			if(normalize){
				V.normalize();
			}
			V.setName(name);
//			V.setHeader(_f.getHeader());
			return V;
		}
	/*
	public CubicVolume solve(CubicVolumeFloat f,float lambda,boolean normalize,int iter){
		int xn=f.getRows();
		int yn=f.getCols();
		int zn=f.getSlices();
		setLabel("Gradient Vector Flow");
		String name=f.getName()+".gvf";
		CubicVolumeFloat V;//=new VectorMatrixX(xn,yn,zn,new Vector3());  
		CubicVolumeFloat Gf;//=new VectorMatrixX(xn,yn,zn,new Vector3());
		CubicVolumeFloat SqrMagf;//=new Matrix3(xn,yn,zn);
		setTotalUnits(5);
		Gf=Gradient3d.doSolve(f);   //compute the gradient of ivol     
		//calculate edge map
		f=Gf.mag();
		f.normalize();
		Gf=Gradient3d.doSolve(f);
		//calculated euclidian norm
		SqrMagf=Gf.mag();
		incrementCompletedUnits();
		V=process(Gf,SqrMagf,lambda,iter);
		if(normalize){
			V.normalize();
		}
		V.setName(name);

		return V;
	}
*/
}
