import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/** P18cyd source-formula reference for default-atlas computeApproxVariances. */
public final class ToadsDefaultApproximateVarianceReferenceHarness {
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String WRAPPER_SHA="d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115";
    private static final String DESCRIPTOR_SHA="cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939";
    private static final int NX=18,NY=18,NZ=18,NC=2,CLASSES=11,APPROX=4;
    private static final float EPSILON=1e-30f;

    private static final class Volume {
        final int nx,ny,nz,nc,datatype;final float[] f;final byte[] b;
        Volume(Path path)throws IOException{
            byte[] raw;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(path));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);raw=out.toByteArray();}
            ByteBuffer h=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);int dims=h.getShort(40);datatype=h.getShort(70);int start=(int)h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dims>3?h.getShort(48):1;
            if(h.getInt(0)!=348||(datatype!=2&&datatype!=16))throw new IOException("unexpected NIfTI datatype");int count=nx*ny*nz*nc;
            if(datatype==16){f=new float[count];b=null;for(int i=0;i<count;i++)f[i]=h.getFloat(start+4*i);}else{f=null;b=new byte[count];for(int i=0;i<count;i++)b[i]=h.get(start+i);}
        }
        int index(int x,int y,int z,int c){return x+nx*(y+ny*(z+nz*c));}
        float value(int x,int y,int z,int c){return f[index(x,y,z,c)];}
        int label(int x,int y,int z,int c){return b[index(x,y,z,c)]&255;}
    }
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static String ints(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static final class Result {float[][] num,den,std;int boundary,duplicates,mismatches,matched;}

    private static Result run(Volume image,Volume field,Volume mem,Volume best,Volume seg,float[][] centers,boolean useField,boolean zeroMembership,boolean removeClass10Support){
        Result r=new Result();r.num=new float[NC][CLASSES];r.den=new float[NC][CLASSES];r.std=new float[NC][CLASSES];
        for(int c=0;c<NC;c++){
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
                int assigned=seg.label(x,y,z,0);if(removeClass10Support&&assigned==10)assigned=9;
                int[] matching=new int[APPROX];int nmatching=0;
                for(int a=0;a<APPROX;a++){
                    int k=best.label(x,y,z,a);if(assigned!=k){r.mismatches++;continue;}r.matched++;matching[nmatching++]=k;
                    boolean boundary=x==0||x==NX-1||y==0||y==NY-1||z==0||z==NZ-1;if(boundary)r.boundary++;
                    float m=mem.value(x,y,z,a);if(zeroMembership&&k==0)m=0.0f;
                    float weight=m*m+EPSILON;
                    float diff;if(useField)diff=field.value(x,y,z,c)*image.value(x,y,z,c)-centers[c][k];else diff=image.value(x,y,z,c)-centers[c][k];
                    r.num[c][k]+=weight*diff*diff;r.den[c][k]+=weight;
                }
                if(nmatching>1)for(int a=1;a<nmatching;a++)if(matching[a]==matching[0]){r.duplicates++;break;}
            }
            for(int k=0;k<CLASSES;k++)if(r.den[c][k]>0.0f)r.std[c][k]=(float)Math.sqrt(r.num[c][k]/r.den[c][k]);
        }
        return r;
    }
    private static int[] flatten(float[][] a){int[] out=new int[a.length*a[0].length];int n=0;for(float[] row:a)for(float v:row)out[n++]=bits(v);return out;}
    public static void main(String[] args)throws Exception{
        Volume image=new Volume(Path.of(args[0])),field=new Volume(Path.of(args[1])),mem=new Volume(Path.of(args[2])),best=new Volume(Path.of(args[3])),seg=new Volume(Path.of(args[4]));float[][] centers=new float[NC][CLASSES];java.util.List<String> rows=Files.readAllLines(Path.of(args[5]));for(int c=0;c<NC;c++){String[] row=rows.get(c).trim().split("\\s+");for(int k=0;k<CLASSES;k++)centers[c][k]=Float.parseFloat(row[k]);}
        Result fielded=run(image,field,mem,best,seg,centers,true,false,false),raw=run(image,field,mem,best,seg,centers,false,false,false),zero=run(image,field,mem,best,seg,centers,true,true,false),empty=run(image,field,mem,best,seg,centers,true,false,true);
        System.out.print("{\"segmentation_sha256\":\""+SEGMENTATION_SHA+"\",\"algorithm_sha256\":\""+ALGORITHM_SHA+"\",\"wrapper_sha256\":\""+WRAPPER_SHA+"\",\"descriptor_sha256\":\""+DESCRIPTOR_SHA+"\",\"shape\":[18,18,18,2],\"classes\":11,\"approx\":4,\"epsilon\":1e-30,\"field_numerator_bits\":"+ints(flatten(fielded.num))+",\"field_denominator_bits\":"+ints(flatten(fielded.den))+",\"field_stddev_bits\":"+ints(flatten(fielded.std))+",\"raw_stddev_bits\":"+ints(flatten(raw.std))+",\"zero_membership_denominator_bits\":"+ints(flatten(zero.den))+",\"zero_membership_stddev_bits\":"+ints(flatten(zero.std))+",\"no_support_class10_stddev_bits\":["+bits(empty.std[0][10])+","+bits(empty.std[1][10])+"] ,\"matching_slots\":"+fielded.matched+",\"mismatching_slots\":"+fielded.mismatches+",\"boundary_matching_slots\":"+fielded.boundary+",\"duplicate_matching_voxels\":"+fielded.duplicates+"}");
    }
}
