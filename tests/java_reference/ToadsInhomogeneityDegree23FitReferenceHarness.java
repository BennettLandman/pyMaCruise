import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

/** Selected-source degree-2/3 fixed polynomial field fit reference. */
public final class ToadsInhomogeneityDegree23FitReferenceHarness {
    private static final String SOURCE_SHA="6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8";
    private static final String JAMA_WRAPPER_SHA="efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526";
    private static final String JAMA_MATRIX_SHA="76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c";
    private static final String JAMA_LU_SHA="9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6";
    private static final String POWER_TABLE_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int STRIDE=3, BACKGROUND=0, NX=18, NY=18, NZ=18;
    private static final class Volume {
        final int nx,ny,nz,nc,datatype;final float[] f;final byte[] b;
        Volume(Path path)throws IOException{
            byte[] raw;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(path));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);raw=out.toByteArray();}
            ByteBuffer h=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);int dims=h.getShort(40);datatype=h.getShort(70);int start=(int)h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dims>3?h.getShort(48):1;
            if(h.getInt(0)!=348||(datatype!=2&&datatype!=16))throw new IOException("unexpected NIfTI datatype");
            int count=nx*ny*nz*nc;if(datatype==16){f=new float[count];b=null;for(int i=0;i<count;i++)f[i]=h.getFloat(start+4*i);}else{f=null;b=new byte[count];for(int i=0;i<count;i++)b[i]=h.get(start+i);}
        }
        int index(int x,int y,int z,int c){return x+nx*(y+ny*(z+nz*c));}
        float value(int x,int y,int z,int c){return f[index(x,y,z,c)];}
        int label(int x,int y,int z,int c){return b[index(x,y,z,c)]&255;}
    }
    private static int bits(float x){return Float.floatToIntBits(x);}
    private static long dbits(double x){return Double.doubleToLongBits(x);}
    private static String ints(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static String longs(double[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(dbits(a[i]));}return s.append(']').toString();}
    private static float[] basis(int x,int y,int z,int degree){
        float rx=x/(float)(NX-1),ry=y/(float)(NY-1),rz=z/(float)(NZ-1);float[] p=new float[degree==2?10:20];
        p[0]=1.0f;p[1]=rx;p[2]=ry;p[3]=rz;
        p[4]=2.0f*rx*rx-1.0f;p[5]=2.0f*ry*ry-1.0f;p[6]=2.0f*rz*rz-1.0f;
        p[7]=p[1]*p[2];p[8]=p[2]*p[3];p[9]=p[3]*p[1];
        if(degree==3){p[10]=4.0f*rx*rx*rx-3.0f*rx;p[11]=4.0f*ry*ry*ry-3.0f*ry;p[12]=4.0f*rz*rz*rz-3.0f*rz;
            p[13]=p[4]*p[2];p[14]=p[1]*p[5];p[15]=p[5]*p[3];p[16]=p[2]*p[6];p[17]=p[6]*p[1];p[18]=p[3]*p[4];p[19]=p[1]*p[2]*p[3];}
        return p;
    }
    private static float optimizedWeight(float image,float centroid){
        float residual=1.0f*image-centroid;float dist=residual*residual;if(dist<=1e-30f)return 1.0f;
        double step=0.000001f;int index=(int)(dist/step+0.5);double cached=((double)index)*step;return (float)(cached/(double)dist);
    }
    private static double[] solveLu(double[][] input,double[] rhs){
        int n=rhs.length;double[][] lu=new double[n][n];for(int i=0;i<n;i++)lu[i]=input[i].clone();int[] piv=new int[n];for(int i=0;i<n;i++)piv[i]=i;
        for(int j=0;j<n;j++){for(int i=0;i<n;i++){double sum=lu[i][j];int kmax=Math.min(i,j);for(int k=0;k<kmax;k++)sum-=lu[i][k]*lu[k][j];lu[i][j]=sum;}
            int p=j;for(int i=j+1;i<n;i++)if(Math.abs(lu[i][j])>Math.abs(lu[p][j]))p=i;if(p!=j){double[] t=lu[p];lu[p]=lu[j];lu[j]=t;int ti=piv[p];piv[p]=piv[j];piv[j]=ti;}
            if(lu[j][j]==0.0)throw new IllegalArgumentException("singular normal matrix");for(int i=j+1;i<n;i++)lu[i][j]/=lu[j][j];}
        double[] x=new double[n];for(int i=0;i<n;i++)x[i]=rhs[piv[i]];for(int k=0;k<n;k++)for(int i=k+1;i<n;i++)x[i]-=x[k]*lu[i][k];
        for(int k=n-1;k>=0;k--){x[k]/=lu[k][k];for(int i=0;i<k;i++)x[i]-=x[k]*lu[i][k];}return x;
    }
    private static String run(Volume image,Volume mem,Volume best,Volume seg,float[] centers,int degree){
        int np=degree==2?10:20;boolean[] mask=new boolean[NX*NY*NZ];double[][] normal=new double[np][np];double[] rhs=new double[np];int maskCount=0,sampleCount=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){boolean active=seg.label(x,y,z,0)!=BACKGROUND;mask[x+NX*(y+NY*z)]=active;if(active)maskCount++;}
        for(int x=0;x<NX;x+=STRIDE)for(int y=0;y<NY;y+=STRIDE)for(int z=0;z<NZ;z+=STRIDE){if(!mask[x+NX*(y+NY*z)])continue;sampleCount++;float[] p=basis(x,y,z,degree);
            for(int a=0;a<mem.nc;a++){int k=best.label(x,y,z,a);if(k==BACKGROUND)continue;float m=mem.value(x,y,z,a);float m2=m*m;float voxel=image.value(x,y,z,0);float w=optimizedWeight(voxel,centers[k]);float vect=m2*w*voxel*centers[k];float matr=m2*w*voxel*voxel;
                for(int n=0;n<np;n++){rhs[n]+=vect*p[n];float rowTerm=matr*p[n];for(int q=0;q<np;q++)normal[n][q]+=(double)(rowTerm*p[q]);}}}
        float volume=(float)(NX*NY*NZ);for(int n=0;n<np;n++){rhs[n]/=volume;for(int m=0;m<np;m++)normal[n][m]/=volume;}
        double[] coeff=solveLu(normal,rhs);float[] field=new float[NX*NY*NZ];Arrays.fill(field,1.0f);float mean=0.0f,den=0.0f;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x+NX*(y+NY*z)]){float[] p=basis(x,y,z,degree);float value=0.0f;for(int n=0;n<np;n++)value=(float)(value+coeff[n]*p[n]);field[x+NX*(y+NY*z)]=value;mean+=value;den++;}
        mean=mean/den;for(int i=0;i<field.length;i++)if(mask[i])field[i]+=1.0f-mean;
        double[] flatNormal=new double[np*np];for(int i=0;i<np;i++)for(int j=0;j<np;j++)flatNormal[i*np+j]=normal[i][j];int[] fieldBits=new int[field.length];for(int i=0;i<field.length;i++)fieldBits[i]=bits(field[i]);
        return "{\"degree\":"+degree+",\"parameter_count\":"+np+",\"mask_count\":"+maskCount+",\"sample_count\":"+sampleCount+",\"rhs_bits\":"+longs(rhs)+",\"normal_matrix_bits\":"+longs(flatNormal)+",\"coefficient_bits\":"+longs(coeff)+",\"field_bits\":"+ints(fieldBits)+",\"mean_before_offset_bits\":"+bits(mean)+"}";
    }
    public static void main(String[] args)throws Exception{
        Volume image=new Volume(Path.of(args[0])),mem=new Volume(Path.of(args[1])),best=new Volume(Path.of(args[2])),seg=new Volume(Path.of(args[3]));float[] centers=new float[11];java.util.List<String> lines=Files.readAllLines(Path.of(args[4]));for(int i=0;i<centers.length;i++)centers[i]=Float.parseFloat(lines.get(i));
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"jama_wrapper_sha256\":\"").append(JAMA_WRAPPER_SHA).append("\",\"jama_matrix_sha256\":\"").append(JAMA_MATRIX_SHA).append("\",\"jama_lu_sha256\":\"").append(JAMA_LU_SHA).append("\",\"power_table_sha256\":\"").append(POWER_TABLE_SHA).append("\",\"shape\":[18,18,18],\"stride\":3,\"background_class\":0,\"degrees\":[");
        out.append(run(image,mem,best,seg,centers,2)).append(',').append(run(image,mem,best,seg,centers,3)).append("]}");System.out.println(out);
    }
}
