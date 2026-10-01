import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.zip.GZIPInputStream;

/** P18cyb synthetic oracle for the selected ratio-gated field-update transition. */
public final class ToadsInhomogeneityTransitionReferenceHarness {
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String CORRECTION_SHA="6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8";
    private static final String JAMA_WRAPPER_SHA="efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526";
    private static final String JAMA_MATRIX_SHA="76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c";
    private static final String JAMA_LU_SHA="9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6";
    private static final String POWER_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int NX=18,NY=18,NZ=18,CHANNELS=2,SLOTS=4,CLASSES=11,BACKGROUND=0;

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
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static long dbits(double v){return Double.doubleToLongBits(v);}
    private static String ints(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static String longs(double[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(dbits(a[i]));}return s.append(']').toString();}
    private static float[] basis(int x,int y,int z,int degree){
        float rx=x/(float)(NX-1),ry=y/(float)(NY-1),rz=z/(float)(NZ-1);float[] p=new float[degree==2?10:20];
        p[0]=1.0f;p[1]=rx;p[2]=ry;p[3]=rz;p[4]=2.0f*rx*rx-1.0f;p[5]=2.0f*ry*ry-1.0f;p[6]=2.0f*rz*rz-1.0f;
        p[7]=p[1]*p[2];p[8]=p[2]*p[3];p[9]=p[3]*p[1];if(degree==3){p[10]=4.0f*rx*rx*rx-3.0f*rx;p[11]=4.0f*ry*ry*ry-3.0f*ry;p[12]=4.0f*rz*rz*rz-3.0f*rz;p[13]=p[4]*p[2];p[14]=p[1]*p[5];p[15]=p[5]*p[3];p[16]=p[2]*p[6];p[17]=p[6]*p[1];p[18]=p[3]*p[4];p[19]=p[1]*p[2]*p[3];}return p;
    }
    private static float optimizedWeight(float field,float image,float centroid){
        float corrected=field*image;float residual=corrected-centroid;float dist=residual*residual;if(dist<=1e-30f)return 1.0f;
        double step=0.000001f;int index=(int)(dist/step+0.5);double cached=((double)index)*step;return (float)(cached/(double)dist);
    }
    private static double[] solveLu(double[][] input,double[] rhs){
        int n=rhs.length;double[][] lu=new double[n][n];for(int i=0;i<n;i++)lu[i]=input[i].clone();int[] piv=new int[n];for(int i=0;i<n;i++)piv[i]=i;
        for(int j=0;j<n;j++){for(int i=0;i<n;i++){double sum=lu[i][j];int kmax=Math.min(i,j);for(int k=0;k<kmax;k++)sum-=lu[i][k]*lu[k][j];lu[i][j]=sum;}int p=j;for(int i=j+1;i<n;i++)if(Math.abs(lu[i][j])>Math.abs(lu[p][j]))p=i;if(p!=j){double[] t=lu[p];lu[p]=lu[j];lu[j]=t;int ti=piv[p];piv[p]=piv[j];piv[j]=ti;}if(lu[j][j]==0.0)throw new IllegalArgumentException("singular normal matrix");for(int i=j+1;i<n;i++)lu[i][j]/=lu[j][j];}
        double[] x=new double[n];for(int i=0;i<n;i++)x[i]=rhs[piv[i]];for(int k=0;k<n;k++)for(int i=k+1;i<n;i++)x[i]-=x[k]*lu[i][k];for(int k=n-1;k>=0;k--){x[k]/=lu[k][k];for(int i=0;i<k;i++)x[i]-=x[k]*lu[i][k];}return x;
    }
    private static String digestField(float[] field)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer bytes=ByteBuffer.allocate(field.length*4).order(ByteOrder.LITTLE_ENDIAN);for(float v:field)bytes.putInt(bits(v));byte[] hash=d.digest(bytes.array());StringBuilder s=new StringBuilder();for(byte b:hash)s.append(String.format("%02x",b&255));return "\""+s+"\"";}
    private static String fit(Volume image,Volume mem,Volume best,Volume seg,float[][] centers,float[] field,int channel,int degree)throws Exception{
        int np=degree==2?10:20;double[][] normal=new double[np][np];double[] rhs=new double[np];int maskCount=0,sampleCount=0;boolean[] mask=new boolean[NX*NY*NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){boolean active=seg.label(x,y,z,0)!=BACKGROUND;mask[x+NX*(y+NY*z)]=active;if(active)maskCount++;}
        for(int x=0;x<NX;x+=3)for(int y=0;y<NY;y+=3)for(int z=0;z<NZ;z+=3){int vi=x+NX*(y+NY*z);if(!mask[vi])continue;sampleCount++;float[] p=basis(x,y,z,degree);float voxel=image.value(x,y,z,channel);
            for(int a=0;a<SLOTS;a++){int k=best.label(x,y,z,a);if(k==BACKGROUND)continue;float m=mem.value(x,y,z,a);float m2=m*m;float w=optimizedWeight(field[vi],voxel,centers[channel][k]);float vect=m2*w*voxel*centers[channel][k];float matr=m2*w*voxel*voxel;
                for(int n=0;n<np;n++){rhs[n]+=vect*p[n];float rowTerm=matr*p[n];for(int q=0;q<np;q++)normal[n][q]+=(double)(rowTerm*p[q]);}}}
        float volume=(float)(NX*NY*NZ);for(int n=0;n<np;n++){rhs[n]/=volume;for(int m=0;m<np;m++)normal[n][m]/=volume;}double[] coeff=solveLu(normal,rhs);float mean=0.0f,den=0.0f;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){int vi=x+NX*(y+NY*z);if(!mask[vi])continue;float[] p=basis(x,y,z,degree);float val=0.0f;for(int n=0;n<np;n++)val=(float)(val+coeff[n]*p[n]);field[vi]=val;mean+=val;den++;}mean=mean/den;for(int i=0;i<field.length;i++)if(mask[i])field[i]+=1.0f-mean;
        double[] flat=new double[np*np];for(int i=0;i<np;i++)for(int j=0;j<np;j++)flat[i*np+j]=normal[i][j];return "{\"degree\":"+degree+",\"mask_count\":"+maskCount+",\"sample_count\":"+sampleCount+",\"rhs_bits\":"+longs(rhs)+",\"normal_matrix_bits\":"+longs(flat)+",\"coefficient_bits\":"+longs(coeff)+",\"field_sha256\":"+digestField(field)+"}";
    }
    private static int next(int degree,int max,float ratio,float low){return degree<max&&ratio<low?degree+1:degree;}
    public static void main(String[] args)throws Exception{
        Volume image=new Volume(Path.of(args[0])),mem=new Volume(Path.of(args[1])),best=new Volume(Path.of(args[2])),seg=new Volume(Path.of(args[3])),initial=new Volume(Path.of(args[4]));float[][] centers=new float[CHANNELS][CLASSES];java.util.List<String> lines=Files.readAllLines(Path.of(args[5]));for(int c=0;c<CHANNELS;c++){String[] row=lines.get(c).trim().split("\\s+");for(int k=0;k<CLASSES;k++)centers[c][k]=Float.parseFloat(row[k]);}
        float[][] fields=new float[CHANNELS][NX*NY*NZ];float[][] attached=new float[CHANNELS][];for(int c=0;c<CHANNELS;c++){for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)fields[c][x+NX*(y+NY*z)]=initial.value(x,y,z,c);attached[c]=fields[c];}
        boolean stableRefs=(attached[0]==fields[0]&&attached[1]==fields[1]);float low=0.5f;int d1=next(1,3,0.4f,low),d2=next(d1,3,0.4f,low);String[] f2=new String[CHANNELS],f3=new String[CHANNELS];for(int c=0;c<CHANNELS;c++)f2[c]=fit(image,mem,best,seg,centers,fields[c],c,d1);for(int c=0;c<CHANNELS;c++)f3[c]=fit(image,mem,best,seg,centers,fields[c],c,d2);
        int[] boundary={next(2,3,Math.nextDown(low),low),next(2,3,low,low),next(2,3,Math.nextUp(low),low),next(3,3,Math.nextDown(low),low)};
        int[] coords={3,4,5,8,7,9,16,14,12};StringBuilder corrected=new StringBuilder("[");for(int i=0;i<coords.length;i+=3){if(i>0)corrected.append(',');int x=coords[i],y=coords[i+1],z=coords[i+2];corrected.append('[');for(int c=0;c<CHANNELS;c++){if(c>0)corrected.append(',');int vi=x+NX*(y+NY*z);corrected.append(bits(fields[c][vi]*image.value(x,y,z,c)));}corrected.append(']');}corrected.append(']');
        StringBuilder out=new StringBuilder("{\"algorithm_sha256\":\"").append(ALGORITHM_SHA).append("\",\"segmentation_sha256\":\"").append(SEGMENTATION_SHA).append("\",\"correction_sha256\":\"").append(CORRECTION_SHA).append("\",\"jama_wrapper_sha256\":\"").append(JAMA_WRAPPER_SHA).append("\",\"jama_matrix_sha256\":\"").append(JAMA_MATRIX_SHA).append("\",\"jama_lu_sha256\":\"").append(JAMA_LU_SHA).append("\",\"power_table_sha256\":\"").append(POWER_SHA).append("\",\"dimensions\":[18,18,18,2],\"approximation_slots\":4,\"low_ratio\":0.5,\"degree_path\":[").append(d1).append(',').append(d2).append("],\"boundary_degrees\":[");for(int i=0;i<boundary.length;i++){if(i>0)out.append(',');out.append(boundary[i]);}out.append("],\"attached_reference_stable\":").append(stableRefs).append(",\"degree2_channel_fits\":[").append(f2[0]).append(',').append(f2[1]).append("],\"degree3_channel_fits\":[").append(f3[0]).append(',').append(f3[1]).append("],\"final_field_sha256\":[").append(digestField(fields[0])).append(',').append(digestField(fields[1])).append("],\"corrected_sample_coordinates\":[[3,4,5],[8,7,9],[16,14,12]],\"corrected_sample_bits\":").append(corrected).append('}');System.out.println(out);
    }
}
