import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/** P18cyc source-formula harness for the selected 11-class computeApproxCentroids branch. */
public final class ToadsDefaultApproximateCentroidReferenceHarness {
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String WRAPPER_SHA="d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115";
    private static final String POWER_VANDY_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final String POWER_COPY_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final String DESCRIPTOR_SHA="cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939";
    private static final int NX=18,NY=18,NZ=18,NC=2,CLASSES=11,APPROX=4;
    private static final float SMOOTH=0.1f,ZERO=1e-10f;
    private static final int[] TYPE={0,1,1,2,2,2,2,2,3,3,3}; // MASK, CSF, GM, WM

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
    private static String ints(int[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);}return s.append(']').toString();}
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static int typeExponent(int type){return type-1;}

    private static Result run(Volume images,Volume field,Volume mems,Volume best,float[][] old,float[][] exponents,boolean emptyClass10){
        float[][] num=new float[NC][CLASSES],den=new float[NC][CLASSES];int zeroFallbacks=0,tableLookups=0,boundarySlots=0,duplicateVoxels=0;
        PowerTable[][] tables=new PowerTable[NC][3];for(int c=0;c<NC;c++)for(int t=0;t<3;t++)tables[c][t]=new PowerTable(0.0f,4.0f,0.000001f,exponents[c][t]);
        for(int c=0;c<NC;c++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            boolean boundary=x==0||x==NX-1||y==0||y==NY-1||z==0||z==NZ-1;int[] seen=new int[APPROX];
            for(int a=0;a<APPROX;a++){
                int k=best.label(x,y,z,a);seen[a]=k;float mem=mems.value(x,y,z,a);if(emptyClass10&&k==10)mem=0.0f;
                if(boundary&&mem!=0.0f)boundarySlots++;
                int type=TYPE[k];float corrected=field.value(x,y,z,c)*images.value(x,y,z,c);float factor=1.0f;
                if(type>0){float dist=(corrected-old[c][k])*(corrected-old[c][k]);if(dist>ZERO){int t=typeExponent(type);factor=(float)(tables[c][t].lookup((double)dist,(double)exponents[c][t])/dist);tableLookups++;}else zeroFallbacks++;}
                float weighted=factor*mem*mem;num[c][k]+=weighted*corrected;den[c][k]+=weighted;
            }
            for(int a=0;a<APPROX;a++)for(int b=a+1;b<APPROX;b++)if(seen[a]==seen[b]){duplicateVoxels++;a=APPROX;break;}
        }
        float[][] updated=new float[NC][CLASSES];for(int c=0;c<NC;c++)for(int k=0;k<CLASSES;k++){updated[c][k]=old[c][k];if(den[c][k]>0.0f)updated[c][k]=(1.0f-SMOOTH)*num[c][k]/den[c][k]+SMOOTH*old[c][k];}
        return new Result(num,den,updated,zeroFallbacks,tableLookups,boundarySlots,duplicateVoxels);
    }
    private static final class Result {
        final float[][] num,den,updated;final int zeroFallbacks,tableLookups,boundarySlots,duplicateVoxels;
        Result(float[][] n,float[][] d,float[][] u,int z,int t,int b,int duplicates){num=n;den=d;updated=u;zeroFallbacks=z;tableLookups=t;boundarySlots=b;duplicateVoxels=duplicates;}
    }
    private static int[] flatten(float[][] a){int[] out=new int[a.length*a[0].length];int n=0;for(float[] row:a)for(float v:row)out[n++]=bits(v);return out;}
    private static int[] flattenFloat(float[][] a){int[] out=new int[a.length*a[0].length];int n=0;for(float[] row:a)for(float v:row)out[n++]=bits(v);return out;}
    public static void main(String[] args)throws Exception{
        Volume images=new Volume(Path.of(args[0])),field=new Volume(Path.of(args[1])),mems=new Volume(Path.of(args[2])),best=new Volume(Path.of(args[3]));float[][] old=new float[NC][CLASSES],exponents=new float[NC][3];java.util.List<String> crows=Files.readAllLines(Path.of(args[4])),erows=Files.readAllLines(Path.of(args[5]));
        for(int c=0;c<NC;c++){String[] a=crows.get(c).trim().split("\\s+"),b=erows.get(c).trim().split("\\s+");for(int k=0;k<CLASSES;k++)old[c][k]=Float.parseFloat(a[k]);for(int t=0;t<3;t++)exponents[c][t]=Float.parseFloat(b[t]);}
        Result result=run(images,field,mems,best,old,exponents,false),empty=run(images,field,mems,best,old,exponents,true);
        System.out.print("{\"segmentation_sha256\":\""+SEGMENTATION_SHA+"\",\"algorithm_sha256\":\""+ALGORITHM_SHA+"\",\"wrapper_sha256\":\""+WRAPPER_SHA+"\",\"power_vandy_sha256\":\""+POWER_VANDY_SHA+"\",\"power_copy_sha256\":\""+POWER_COPY_SHA+"\",\"descriptor_sha256\":\""+DESCRIPTOR_SHA+"\",\"shape\":[18,18,18,2],\"classes\":11,\"approx\":4,\"centroid_mode\":\"damped\",\"smoothness\":0.1,\"zero\":1e-10,\"numerator_bits\":"+ints(flatten(result.num))+",\"denominator_bits\":"+ints(flatten(result.den))+",\"centroid_bits\":"+ints(flatten(result.updated))+",\"zero_fallbacks\":"+result.zeroFallbacks+",\"table_lookups\":"+result.tableLookups+",\"boundary_slots\":"+result.boundarySlots+",\"duplicate_candidate_voxels\":"+result.duplicateVoxels+",\"empty_class10_denominator_bits\":"+bits(empty.den[0][10])+",\"empty_class10_centroid_bits\":"+bits(empty.updated[0][10])+",\"empty_class10_centroid_bits_ch1\":"+bits(empty.updated[1][10])+"}");
    }
}
