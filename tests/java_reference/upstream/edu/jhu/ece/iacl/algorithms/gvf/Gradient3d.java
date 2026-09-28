package edu.jhu.ece.iacl.algorithms.gvf;
import javax.vecmath.Vector3f;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.MaskVolume18;
/**
 * Compute gradient for volume using differentiator df[i]=0.5*(f[i+1]-f[i-1])
 * @author Blake Lucas
 *
 */
public class Gradient3d {
	public static ImageDataFloat doSolve(ImageDataFloat A) {
		int xn = A.getRows();
		int yn = A.getCols();
		int zn = A.getSlices();
		int x, y, z;
		ImageDataFloat B = new ImageDataFloat(xn, yn, zn, 3);
		float[][][][] Bmat=B.toArray4d();
		float[][][] Amat=A.toArray3d();
		/* compute the x component gradient */
		for (y = 0; y < yn; y++) {
			for (z = 0; z < zn; z++) {
				/* deal with boundary conditions */
				Bmat[0][y][z][0]=( Amat[1][y][z] - Amat[0][y][z]);
				Bmat[xn - 1][y][z][0]=( Amat[xn - 1][y][z]
						- Amat[xn - 2][y][z]);
				for (x = 1; x < xn - 1; x++) {
					Bmat[x][y][z][0]=( (Amat[x + 1][y][z] - Amat[x - 1][y][z]) *0.5f);
				}
			}
		}

		/* compute the y component gradient */
		for (x = 0; x < xn; x++) {
			for (z = 0; z < zn; z++) {
				/* deal with boundary conditions */
				Bmat[x][0][z][1]=( Amat[x][1][z] - Amat[x][0][z]);
				Bmat[x][yn - 1][z][1]=( Amat[x][yn - 1][z]
						- Amat[x][yn - 2][z]);

				for (y = 1; y < yn - 1; y++) {
					Bmat[x][y][z][1]=( (Amat[x][y + 1][z] - Amat[x][y - 1][z]) *0.5f);
				}
			}
		}
		/* compute the z component gradient */
		for (x = 0; x < xn; x++) {
			for (y = 0; y < yn; y++) {
				/* deal with boundary conditions */
				Bmat[x][y][0][2]=( Amat[x][y][1] - Amat[x][y][0]);
				Bmat[x][y][zn - 1][2]=( Amat[x][y][zn - 1]
						- Amat[x][y][zn - 2]);

				for (z = 1; z < zn - 1; z++) {
					Bmat[x][y][z][2]=((Amat[x][y][z + 1] - Amat[x][y][z - 1]) *0.5f);
				}
			}
		}
		return B;
	}
	public static ImageDataFloat doSolve18(ImageDataFloat A) {
		int xn = A.getRows();
		int yn = A.getCols();
		int zn = A.getSlices();
		int x, y, z;
		ImageDataFloat B = new ImageDataFloat(xn, yn, zn, 3);
		float[][][][] Bmat=B.toArray4d();
		float[][][] Amat=A.toArray3d();
		/* compute the z component gradient */
		final float h1=1.0f/(2*(float)Math.sqrt(2));
		final float h2=1.0f/(2*(float)Math.sqrt(2));
		final float h3=1.0f/(2*(float)Math.sqrt(2));
		final float h4=1.0f/(2*(float)Math.sqrt(2));
		MaskVolume18 mask=new MaskVolume18();
		Vector3f[] vec=new Vector3f[18];
		byte[] xnbhd=mask.getNeighborsX();
		byte[] ynbhd=mask.getNeighborsY();
		byte[] znbhd=mask.getNeighborsZ();
		for(int i=0;i<18;i++){
			vec[i]=new Vector3f((float)xnbhd[i],(float)ynbhd[i],(float)znbhd[i]);
			vec[i].normalize();
		}
		Vector3f v;
		for (x = 1; x < xn-1; x++) {
			for (y = 1; y < yn-1; y++) {
				for(z=1; z < zn-1; z++){
					for(int i=0;i<18;i++){
						v=vec[i];
						Bmat[x][y][z][0]+=0.5*v.x*Amat[x+xnbhd[i]][y+ynbhd[i]][z+znbhd[i]];
						Bmat[x][y][z][1]+=0.5*v.y*Amat[x+xnbhd[i]][y+ynbhd[i]][z+znbhd[i]];
						Bmat[x][y][z][2]+=0.5*v.z*Amat[x+xnbhd[i]][y+ynbhd[i]][z+znbhd[i]];
					}
				}
			}
		}

		return B;
	}
	public static ImageDataFloat doSolve(ImageDataFloat A,byte[][][] mask){
		double dx=0,dy=0,dz=0;
		double d1=0,d2=0;
		double uxx=0,uyy=0,uzz=0;
		int xn = A.getRows();
		int yn = A.getCols();
		int zn = A.getSlices();
		int x, y, z;
		ImageDataFloat B = new ImageDataFloat(xn, yn, zn, 3);
		float[][][][] Bmat=B.toArray4d();
		float[][][] Amat=A.toArray3d();
		for(int i=0;i<xn;i++){
			for(int j=0;j<yn;j++){
				for(int k=0;k<zn;k++){
					dx=dy=dz=0;
					uxx=uyy=uzz=0;
					if(i+1<xn&&mask[i+1][j][k]==1){
						uxx=Amat[i+1][j][k]-Amat[i][j][k];
						dx++;
					} 
					if(i>0&&mask[i-1][j][k]==1){
						uxx+=Amat[i][j][k]-Amat[i-1][j][k];
						dx++;
					}
					if(j+1<yn&&mask[i][j+1][k]==1){
						uyy=Amat[i][j+1][k]-Amat[i][j][k];
						dy++;
					} 
					if(j>0&&mask[i][j-1][k]==1){
						uyy+=Amat[i][j][k]-Amat[i][j-1][k];
						dy++;
					}
					if(k+1<zn&&mask[i][j][k+1]==1){
						uzz=Amat[i][j][k+1]-Amat[i][j][k];
						dz++;
					} 
					if(k>0&&mask[i][j][k-1]==1){
						uzz+=Amat[i][j][k]-Amat[i][j][k-1];
						dz++;
					}
					dx=Math.max(dx, 1);
					dy=Math.max(dy, 1);
					dz=Math.max(dz, 1);
					Bmat[i][j][k][0]=(float)(uxx/dx);
					Bmat[i][j][k][1]=(float)(uyy/dy);
					Bmat[i][j][k][2]=(float)(uzz/dz);
				}
			}
		}
		return B;
	}
}
