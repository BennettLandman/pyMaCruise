import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected atlas-aware initialApproxMemberships excerpt for default TOADS. */
public final class ToadsDefaultApproximateMembershipReferenceHarness {
    private static final String SOURCE_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String DESCRIPTOR_SHA="cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939";
    private static final int APPROX=3;
    private static final class Volume {
        final int nx,ny,nz,nc; final float[][][][] data;
        Volume(Path file)throws IOException{
            byte[] bytes;
            try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(file));ByteArrayOutputStream out=new ByteArrayOutputStream()){
                in.transferTo(out);bytes=out.toByteArray();
            }
            ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
            int dims=h.getShort(40),datatype=h.getShort(70);float offset=h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dims>3?h.getShort(48):1;
            if(h.getInt(0)!=348||datatype!=16||nx<1||ny<1||nz<1||nc<1)throw new IOException("expected float32 NIfTI-1");
            data=new float[nx][ny][nz][nc];int start=(int)offset;
            for(int c=0;c<nc;c++)for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
                int i=x+nx*(y+ny*(z+nz*c));data[x][y][z][c]=h.getFloat(start+4*i);
            }
        }
    }
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String floats(float[][][][] values){
        StringBuilder out=new StringBuilder("[");boolean first=true;
        for(float[][][] x:values)for(float[][] y:x)for(float[] z:y)for(float v:z){if(!first)out.append(',');first=false;out.append(bits(v));}
        return out.append(']').toString();
    }
    private static String bytes(byte[][][][] values){
        StringBuilder out=new StringBuilder("[");boolean first=true;
        for(byte[][][] x:values)for(byte[][] y:x)for(byte[] z:y)for(byte v:z){if(!first)out.append(',');first=false;out.append(v&255);}
        return out.append(']').toString();
    }
    private static String strings(List<String> values){
        StringBuilder out=new StringBuilder("[");
        for(int i=0;i<values.size();i++){if(i>0)out.append(',');out.append('"').append(values.get(i)).append('"');}
        return out.append(']').toString();
    }
    private static List<String> classTypes(Path descriptor)throws IOException{
        List<String> lines=Files.readAllLines(descriptor);int head=-1;
        for(int i=0;i<lines.size();i++)if(lines.get(i).startsWith("Structures:")){head=i;break;}
        if(head<0)throw new IOException("descriptor lacks Structures header");
        int count=Integer.parseInt(lines.get(head).split(":",2)[1].trim());List<String> types=new ArrayList<>();
        for(int k=0;k<count;k++){
            String[] fields=lines.get(head+1+k).trim().split("\\s+");String tag=fields[fields.length-1];
            types.add(tag.equals("mask")?"MASK":tag.toUpperCase());
        }
        return types;
    }
    private static float[][] centroids(int channels,int classes){
        float[][] out=new float[channels][classes];
        for(int c=0;c<channels;c++)for(int k=0;k<classes;k++)out[c][k]=.19f+.021f*k+.013f*c;
        return out;
    }
    private static void selectedInitialApprox(Volume image,Volume shape,float[][] center,float[][][][] memberships,byte[][][][] best){
        int classes=shape.nc,channels=image.nc;
        for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
            float[] num=new float[classes];float den=0.0f;
            for(int k=0;k<classes;k++){
                num[k]=0.0f;
                for(int c=0;c<channels;c++){
                    float delta=image.data[x][y][z][c]-center[c][k];num[k]+=delta*delta;
                }
                for(int m=0;m<classes;m++)if(m!=k){float prior=shape.data[x][y][z][m];num[k]+=prior*prior;}
                if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
                den+=num[k];
            }
            if(den>0.0f)for(int k=0;k<classes;k++)num[k]/=den;
            float[] localShape=new float[classes];
            for(int k=0;k<classes;k++)localShape[k]=shape.data[x][y][z][k];
            Numerics.bestIndex(best[x][y][z],memberships[x][y][z],num,APPROX,Numerics.bestIndex(localShape));
        }
    }
    public static void main(String[] args)throws Exception{
        List<String> types=classTypes(Path.of(args[2]));
        List<String> expected=Arrays.asList("MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_WM","OPTIMIZED_WM","OPTIMIZED_WM");
        if(!types.equals(expected))throw new IllegalArgumentException("default TOADS descriptor order mismatch: "+types);
        Volume image=new Volume(Path.of(args[0])),shape=new Volume(Path.of(args[1]));
        if(image.nx!=shape.nx||image.ny!=shape.ny||image.nz!=shape.nz||image.nc!=2||shape.nc!=11)throw new IllegalArgumentException("synthetic fixture grid/class mismatch");
        float[][][][] memberships=new float[image.nx][image.ny][image.nz][APPROX];
        byte[][][][] best=new byte[image.nx][image.ny][image.nz][APPROX];
        selectedInitialApprox(image,shape,centroids(image.nc,shape.nc),memberships,best);
        System.out.println("{\"source_sha256\":\""+SOURCE_SHA+"\",\"descriptor_sha256\":\""+DESCRIPTOR_SHA+
            "\",\"descriptor_class_types\":"+strings(types)+",\"approx\":"+APPROX+
            ",\"initial_memberships_bits\":"+floats(memberships)+",\"best_class_indices\":"+bytes(best)+"}");
    }
}
