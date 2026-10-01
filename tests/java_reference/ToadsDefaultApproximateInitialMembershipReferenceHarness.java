import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected no-atlas initialApproxMemberships excerpt for default TOADS. */
public final class ToadsDefaultApproximateInitialMembershipReferenceHarness {
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String LESION_ALGORITHM_SHA="556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final String LESION_PREPROCESS_SHA="72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08";
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
    private static Object[] selectedVoxel(float[] image,float[][] centroid){
        int classes=centroid[0].length;float den=0.0f;float[] num=new float[classes],grouping=new float[classes];
        for(int k=0;k<classes;k++){
            num[k]=0.0f;
            for(int c=0;c<image.length;c++){float delta=image[c]-centroid[c][k];num[k]+=delta*delta;}
            if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
            den+=num[k];grouping[k]=0.0f;
            for(int m=0;m<classes;m++)if(centroid[0][m]==centroid[0][k])grouping[k]++;
        }
        if(den>0.0f)for(int k=0;k<classes;k++)num[k]=grouping[k]*num[k]/den;
        byte[] best=new byte[APPROX];float[] memberships=new float[APPROX];
        Numerics.bestIndex(best,memberships,num,APPROX);
        return new Object[]{best,memberships};
    }
    private static String vector(float[] values){StringBuilder out=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.floatToIntBits(values[i]));}return out.append(']').toString();}
    private static String byteVector(byte[] values){StringBuilder out=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]&255);}return out.append(']').toString();}
    private static String sha256(Path file)throws Exception{
        byte[] digest=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(file));StringBuilder out=new StringBuilder();
        for(byte value:digest)out.append(String.format("%02x",value&255));return out.toString();
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
        StringBuilder output=new StringBuilder("{\"algorithm_source_sha256\":\"").append(ALGORITHM_SHA).append("\",\"preprocess_source_sha256\":\"").append(PREPROCESS_SHA)
            .append("\",\"image_sha256\":\"").append(sha256(Path.of(args[0])))
            .append("\",\"descriptor_sha256\":\"").append(DESCRIPTOR_SHA).append("\",\"lesion_algorithm_source_sha256\":\"").append(LESION_ALGORITHM_SHA)
            .append("\",\"lesion_preprocess_source_sha256\":\"").append(LESION_PREPROCESS_SHA)
            .append("\",\"selected_wrapper\":\"AlgorithmMultipleToads\",\"selected_overload\":\"initialMemberships(memberships,best,approx,centroids,classes)\"")
            .append(",\"wrapper_call_order\":[\"AlgorithmMultipleToads.runAlgorithm\",\"preprocess.initialMemberships(memberships,best,approx,centroids,classes)\",\"atlas.setImageInfo\",\"atlas.initShapeRegistration\",\"atlas.updateShapeRegistration\",\"atlas.registerShapesPyramid\"]")
            .append(",\"dedicated_lesion_wrapper_call\":\"AlgorithmLesionToads -> LesionToadPreprocess.initialMemberships(centroids,classes): separate full-class method, no approx parameter\"")
            .append(",\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],\"approx\":"+APPROX)
            .append(",\"has_atlas_shape_penalty\":false,\"has_forced_candidate_slot\":false")
            .append(",\"initial_memberships_bits\":").append(bits(mem)).append(",\"best_class_indices\":").append(indexes(best));
        Object[] zero=selectedVoxel(new float[]{0.0f},new float[1][11]);
        output.append(",\"zero_cost_tie\":{\"indices\":").append(byteVector((byte[])zero[0])).append(",\"membership_bits\":").append(vector((float[])zero[1])).append('}');
        Object[] near=selectedVoxel(new float[]{1.0e-11f},new float[1][11]);
        output.append(",\"below_cutoff_tie\":{\"indices\":").append(byteVector((byte[])near[0])).append(",\"membership_bits\":").append(vector((float[])near[1])).append('}');
        float[][] multiplicityCentroids=new float[1][11];multiplicityCentroids[0][3]=1.0e-10f;
        for(int k=4;k<11;k++)multiplicityCentroids[0][k]=2.0e-10f;
        Object[] multiplicity=selectedVoxel(new float[]{0.0f},multiplicityCentroids);
        output.append(",\"duplicate_centroid_multiplicity\":{\"indices\":").append(byteVector((byte[])multiplicity[0])).append(",\"membership_bits\":").append(vector((float[])multiplicity[1])).append('}').append('}');
        System.out.println(output);
    }
}
