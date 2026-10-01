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

/** Selected no-atlas initialApproxMemberships excerpt for default TOADS. */
public final class ToadsDefaultApproximateInitialMembershipReferenceHarness {
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
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
    private static float[][] centroids(int channels,int classes){
        float[][] out=new float[channels][classes];
        for(int c=0;c<channels;c++)for(int k=0;k<classes;k++)out[c][k]=.19f+.021f*k+.013f*c;
        out[0][2]=out[0][1];
        return out;
    }
    private static void initialApproxMemberships(Volume image,float[][] centroid,float[][][][] mem,byte[][][][] best){
        int classes=centroid[0].length;
        for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
            float den=0.0f;float[] num=new float[classes],grouping=new float[classes];
            for(int k=0;k<classes;k++){
                num[k]=0.0f;
                for(int c=0;c<image.nc;c++){
                    float delta=image.data[x][y][z][c]-centroid[c][k];num[k]+=delta*delta;
                }
                if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
                den+=num[k];
                grouping[k]=0.0f;
                for(int m=0;m<classes;m++)if(centroid[0][m]==centroid[0][k])grouping[k]++;
            }
            if(den>0.0f)for(int k=0;k<classes;k++)num[k]=grouping[k]*num[k]/den;
            Numerics.bestIndex(best[x][y][z],mem[x][y][z],num,APPROX);
        }
    }
    private static String bits(float[][][][] a){
        StringBuilder out=new StringBuilder("[");boolean first=true;
        for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){if(!first)out.append(',');first=false;out.append(Float.floatToIntBits(v));}
        return out.append(']').toString();
    }
    private static String indexes(byte[][][][] a){
        StringBuilder out=new StringBuilder("[");boolean first=true;
        for(byte[][][] x:a)for(byte[][] y:x)for(byte[] z:y)for(byte v:z){if(!first)out.append(',');first=false;out.append(v&255);}
        return out.append(']').toString();
    }
    private static List<String> classTypes(Path descriptor)throws IOException{
        List<String> lines=Files.readAllLines(descriptor);int head=-1;
        for(int i=0;i<lines.size();i++)if(lines.get(i).startsWith("Structures:")){head=i;break;}
        if(head<0)throw new IOException("descriptor lacks Structures header");
        int count=Integer.parseInt(lines.get(head).split(":",2)[1].trim());List<String> out=new ArrayList<>();
        for(int i=0;i<count;i++){
            String[] fields=lines.get(head+1+i).trim().split("\\s+");String tag=fields[fields.length-1];
            if(tag.equals("mask"))out.add("MASK");else out.add(tag.toUpperCase());
        }
        return out;
    }
    public static void main(String[] args)throws Exception{
        List<String> types=classTypes(Path.of(args[1]));
        List<String> expected=Arrays.asList("MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_WM","OPTIMIZED_WM","OPTIMIZED_WM");
        if(!types.equals(expected))throw new IllegalArgumentException("unexpected default TOADS class order: "+types);
        Volume image=new Volume(Path.of(args[0]));if(image.nc!=2)throw new IllegalArgumentException("expected two-channel phantom");
        float[][][][] mem=new float[image.nx][image.ny][image.nz][APPROX];byte[][][][] best=new byte[image.nx][image.ny][image.nz][APPROX];
        initialApproxMemberships(image,centroids(image.nc,types.size()),mem,best);
        System.out.println("{\"algorithm_source_sha256\":\""+ALGORITHM_SHA+"\",\"preprocess_source_sha256\":\""+PREPROCESS_SHA+
            "\",\"descriptor_sha256\":\""+DESCRIPTOR_SHA+"\",\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],\"approx\":"+APPROX+
            ",\"initial_memberships_bits\":"+bits(mem)+",\"best_class_indices\":"+indexes(best)+'}');
    }
}
