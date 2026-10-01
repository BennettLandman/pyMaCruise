import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/**
 * Selected-source formula harness for InhomogeneityCorrection2's degree-one
 * approximate polynomial field. The 4x4 solve reproduces the square Jama LU
 * path used by Jama.Matrix.solve; all data and field accumulation are float32.
 */
public final class ToadsInhomogeneityDegreeOneFitReferenceHarness {
    private static final String SOURCE_SHA="6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8";
    private static final String JAMA_WRAPPER_SHA="efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526";
    private static final String JAMA_MATRIX_SHA="76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c";
    private static final String JAMA_LU_SHA="9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6";
    private static final String POWER_TABLE_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int STRIDE=3, DEGREE=1, NP=4, BACKGROUND=0;
    private static final class Volume {
        final int nx,ny,nz,nc,datatype;final float[] f;final byte[] b;
        Volume(Path path)throws IOException{
            byte[] bytes;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(path));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);bytes=out.toByteArray();}
            ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);int dimensions=h.getShort(40);datatype=h.getShort(70);float offset=h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dimensions>3?h.getShort(48):1;int count=nx*ny*nz*nc,start=(int)offset;
            if(h.getInt(0)!=348||(datatype!=2&&datatype!=16))throw new IOException("expected float32 or uint8 NIfTI-1");
            if(datatype==16){f=new float[count];b=null;for(int i=0;i<count;i++)f[i]=h.getFloat(start+4*i);}else{f=null;b=new byte[count];for(int i=0;i<count;i++)b[i]=h.get(start+i);}
        }
        int index(int x,int y,int z,int c){return x+nx*(y+ny*(z+nz*c));}
        float value(int x,int y,int z,int c){return f[index(x,y,z,c)];}
        int label(int x,int y,int z,int c){return b[index(x,y,z,c)]&255;}
    }
    private static int bits(float x){return Float.floatToIntBits(x);}
    private static long dbits(double x){return Double.doubleToLongBits(x);}
    private static String intArray(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static String longArray(double[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(dbits(a[i]));}return s.append(']').toString();}
    private static String floatArray(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static float[] basis(int x,int y,int z,int nx,int ny,int nz){return new float[]{1.0f,x/(float)(nx-1),y/(float)(ny-1),z/(float)(nz-1)};}
    private static float optimizedWeight(float image,float centroid){
        float residual=1.0f*image-centroid;float dist=residual*residual;
        if(dist<=1e-30f)return 1.0f;
        double step=0.000001f;int index=(int)(dist/step+0.5);double cached=((double)index)*step;
        return (float)(cached/(double)dist);
    }
    private static double[] solveLu(double[][] input,double[] rhs){
        int n=rhs.length;double[][] lu=new double[n][n];for(int i=0;i<n;i++)lu[i]=input[i].clone();int[] piv=new int[n];for(int i=0;i<n;i++)piv[i]=i;
        for(int j=0;j<n;j++){
            for(int i=0;i<n;i++){double sum=lu[i][j];int kmax=Math.min(i,j);for(int k=0;k<kmax;k++)sum-=lu[i][k]*lu[k][j];lu[i][j]=sum;}
            int p=j;for(int i=j+1;i<n;i++)if(Math.abs(lu[i][j])>Math.abs(lu[p][j]))p=i;
            if(p!=j){double[] t=lu[p];lu[p]=lu[j];lu[j]=t;int ti=piv[p];piv[p]=piv[j];piv[j]=ti;}
            if(j<n&&lu[j][j]!=0.0)for(int i=j+1;i<n;i++)lu[i][j]/=lu[j][j];
        }
        double[] x=new double[n];for(int i=0;i<n;i++)x[i]=rhs[piv[i]];
        for(int k=0;k<n;k++)for(int i=k+1;i<n;i++)x[i]-=x[k]*lu[i][k];
        for(int k=n-1;k>=0;k--){x[k]/=lu[k][k];for(int i=0;i<k;i++)x[i]-=x[k]*lu[i][k];}
        return x;
    }
    private static String run(Volume image,Volume mem,Volume best,Volume seg,float[] centers){
        int nx=image.nx,ny=image.ny,nz=image.nz;boolean[] mask=new boolean[nx*ny*nz];double[][] factor=new double[NP][NP];double[] rhs=new double[NP];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)mask[x+nx*(y+ny*z)]=seg.label(x,y,z,0)!=BACKGROUND;
        for(int x=0;x<nx;x+=STRIDE)for(int y=0;y<ny;y+=STRIDE)for(int z=0;z<nz;z+=STRIDE){
            if(!mask[x+nx*(y+ny*z)])continue;float[] pol=basis(x,y,z,nx,ny,nz);for(int a=0;a<mem.nc;a++){
                int k=best.label(x,y,z,a);if(k==BACKGROUND)continue;float m=mem.value(x,y,z,a);float m2=m*m;float voxel=image.value(x,y,z,0);float weight=optimizedWeight(voxel,centers[k]);
                float vector=m2*weight*voxel*centers[k];float matrix=m2*weight*voxel*voxel;
                for(int n=0;n<NP;n++){rhs[n]+=vector*pol[n];for(int q=0;q<NP;q++)factor[n][q]+=matrix*pol[n]*pol[q];}
            }
        }
        float volume=(float)(nx*ny*nz);for(int n=0;n<NP;n++){rhs[n]/=volume;for(int m=0;m<NP;m++)factor[n][m]/=volume;}
        double[] coeff=solveLu(factor,rhs);float[] field=new float[nx*ny*nz];java.util.Arrays.fill(field,1.0f);float mean=0.0f,den=0.0f;
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)if(mask[x+nx*(y+ny*z)]){
            float[] pol=basis(x,y,z,nx,ny,nz);float value=0.0f;for(int n=0;n<NP;n++)value+=(float)(coeff[n]*pol[n]);field[x+nx*(y+ny*z)]=value;mean+=value;den++;
        }
        mean=mean/den;for(int i=0;i<field.length;i++)if(mask[i])field[i]+=1.0f-mean;
        double[] flatFactor=new double[NP*NP];for(int i=0;i<NP;i++)for(int j=0;j<NP;j++)flatFactor[i*NP+j]=factor[i][j];
        int nmask=0;for(boolean v:mask)if(v)nmask++;
        return "{\"source_sha256\":\""+SOURCE_SHA+"\",\"jama_wrapper_sha256\":\""+JAMA_WRAPPER_SHA+"\",\"jama_matrix_sha256\":\""+JAMA_MATRIX_SHA+"\",\"jama_lu_sha256\":\""+JAMA_LU_SHA+"\",\"power_table_sha256\":\""+POWER_TABLE_SHA+"\",\"degree\":"+DEGREE+",\"stride\":"+STRIDE+",\"mask_count\":"+nmask+",\"rhs_bits\":"+longArray(rhs)+",\"normal_matrix_bits\":"+longArray(flatFactor)+",\"coefficient_bits\":"+longArray(coeff)+",\"field_bits\":"+floatArray(field)+"}";
    }
    public static void main(String[] args)throws Exception{
        float[] centers=new float[11];for(int i=0;i<centers.length;i++)centers[i]=Float.parseFloat(Files.readAllLines(Path.of(args[4])).get(i));
        System.out.println(run(new Volume(Path.of(args[0])),new Volume(Path.of(args[1])),new Volume(Path.of(args[2])),new Volume(Path.of(args[3])),centers));
    }
}
