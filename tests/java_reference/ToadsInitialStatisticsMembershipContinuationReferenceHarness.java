import edu.jhmi.rad.medic.utilities.PowerTable;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.Arrays;

/**
 * Fixed-input P18cyn Java continuation: post-thinning exact centroids,
 * assigned-class standard deviations, centroid similarity, then memberships.
 * It reuses the independently pinned P18cyn source harness without changing
 * that harness or its existing source hash.
 */
public final class ToadsInitialStatisticsMembershipContinuationReferenceHarness {
    private static final int C=11;
    private static Object get(Object target,String name)throws Exception{
        Field f=target.getClass().getDeclaredField(name);f.setAccessible(true);return f.get(target);
    }
    private static void set(Object target,String name,Object value)throws Exception{
        Field f=target.getClass().getDeclaredField(name);f.setAccessible(true);f.set(target,value);
    }
    private static Object getStatic(String name)throws Exception{
        Field f=ToadsInitialThinningStartupReferenceHarness.class.getDeclaredField(name);f.setAccessible(true);return f.get(null);
    }
    private static String hex(byte[] bytes){StringBuilder s=new StringBuilder();for(byte b:bytes)s.append(String.format("%02x",b&255));return s.toString();}
    private static void put(MessageDigest d,float value){d.update(ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN).putInt(Float.floatToIntBits(value)).array());}
    private static String hash(float[][][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(float[][][] p:a)for(float[][] r:p)for(float[] v:r)for(float q:v)put(d,q);return hex(d.digest());}
    private static String hash(float[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(float[][] p:a)for(float[] r:p)for(float q:r)put(d,q);return hex(d.digest());}
    private static String hash(float[][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(float[] r:a)for(float q:r)put(d,q);return hex(d.digest());}
    private static String hash(byte[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(byte[][] p:a)for(byte[] r:p)d.update(r);return hex(d.digest());}
    private static String bits(float[][] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append('[');for(int j=0;j<a[i].length;j++){if(j>0)s.append(',');s.append(Float.floatToIntBits(a[i][j]));}s.append(']');}return s.append(']').toString();}

    private static void refreshStatistics(Object h)throws Exception{
        float[][][][] images=(float[][][][])get(h,"image"), fields=(float[][][][])get(h,"attachedFieldByChannel");
        float[][][][] memberships=(float[][][][])get(h,"memberships");
        byte[][][] segmentation=(byte[][][])get(h,"segmentation");
        float[][] centers=(float[][])getStatic("CENTROID"), deviations=(float[][])getStatic("VARIANCE");
        float[][] modalityExponent=new float[][]{{.78f,.91f},{1.16f,1.29f},{1.55f,1.73f}};
        float[][] old=new float[2][C];
        for(int c=0;c<2;c++)System.arraycopy(centers[c],0,old[c],0,C);
        for(int c=0;c<2;c++){
            PowerTable[] powers={new PowerTable(0,4,.000001f,modalityExponent[0][c]),new PowerTable(0,4,.000001f,modalityExponent[1][c]),new PowerTable(0,4,.000001f,modalityExponent[2][c])};
            for(int k=0;k<C;k++){
                float numerator=0.0f,denominator=0.0f;
                for(int x=1;x<images.length-1;x++)for(int y=1;y<images[x].length-1;y++)for(int z=1;z<images[x][y].length-1;z++){
                    float corrected=fields[c][x][y][z]*images[x][y][z][c];
                    float mem=memberships[x][y][z][k], weight;
                    if(k==0)weight=mem*mem;
                    else{
                        float residual=corrected-old[c][k], distance=residual*residual, factor=1.0f;
                        if(distance>1.0e-10f)factor=(float)(powers[k<=2?0:k<=7?1:2].lookup((double)distance,(double)modalityExponent[k<=2?0:k<=7?1:2][c])/(double)distance);
                        weight=(factor*mem)*mem;
                    }
                    numerator+=weight*corrected;denominator+=weight;
                }
                if(denominator>0.0f)centers[c][k]=((1.0f-.1f)*numerator)/denominator+.1f*old[c][k];
            }
        }
        for(int c=0;c<2;c++)for(int k=0;k<C;k++){
            float numerator=0.0f,denominator=0.0f;
            for(int x=0;x<images.length;x++)for(int y=0;y<images[x].length;y++)for(int z=0;z<images[x][y].length;z++)if((segmentation[x][y][z]&255)==k){
                float mem=memberships[x][y][z][k],weight=mem*mem+1e-30f;
                float corrected=fields[c][x][y][z]*images[x][y][z][c],delta=corrected-centers[c][k];
                numerator+=(weight*delta)*delta;denominator+=weight;
            }
            deviations[c][k]=denominator>0.0f?(float)Math.sqrt(numerator/denominator):0.0f;
        }
        float[][] similarity=(float[][])getStatic("SIM");
        float[] imageVariance=(float[])getStatic("IMAGE_VARIANCE");
        float scaleFactor=(1.0f-.05f)/.05f;
        for(int k=0;k<C;k++)for(int m=0;m<C;m++){
            if(k==m){similarity[k][m]=0.0f;continue;}
            float distance=0.0f;
            for(int c=0;c<2;c++){
                float delta=centers[c][k]-centers[c][m],sum=deviations[c][k]+deviations[c][m];
                float term=imageVariance[c]*delta;term*=delta;term/=sum;term/=sum;distance+=term;
            }
            similarity[k][m]=1.0f/(1.0f+scaleFactor*distance);
        }
    }

    public static void main(String[] args)throws Exception{
        if(args.length!=1)throw new IllegalArgumentException("fixture directory required");
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        Constructor<?> ctor=ToadsInitialThinningStartupReferenceHarness.class.getDeclaredConstructor(String.class);ctor.setAccessible(true);
        Object h=ctor.newInstance(args[0]);
        refreshStatistics(h);
        float[][][][] current=(float[][][][])get(h,"memberships");
        String beforeMembership=hash(current);
        set(h,"initial",current);
        set(h,"objective",0.0f);
        Method next=ToadsInitialThinningStartupReferenceHarness.class.getDeclaredMethod("computeMemberships");next.setAccessible(true);next.invoke(h);
        float[][][][] energy=(float[][][][])get(h,"energy"), memberships=(float[][][][])get(h,"memberships");
        float[][][] denominator=(float[][][])get(h,"denominator");
        float[][] centers=(float[][])getStatic("CENTROID"),deviations=(float[][])getStatic("VARIANCE"),similarity=(float[][])getStatic("SIM");
        byte[][][] segmentation=(byte[][][])get(h,"segmentation");
        float objective=(float)get(h,"objective");
        System.out.print("{\"source\":\"P18cyn startup/first-thinning -> computeExactCentroids -> computeExactVariances -> centroid similarity -> next computeExactMemberships\"");
        System.out.print(",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\"");
        System.out.print(",\"source_order\":[\"P18cyn:startup membership and first thinning\",\"computeExactCentroids\",\"computeExactVariances\",\"computeCentroidSimilarity\",\"computeExactMemberships\"]");
        System.out.print(",\"startup_membership_sha256\":\""+beforeMembership+"\",\"post_thinning_segmentation_sha256\":\""+hash(segmentation)+"\"");
        System.out.print(",\"updated_centroids_bits\":"+bits(centers)+",\"standard_deviation_bits\":"+bits(deviations));
        System.out.print(",\"centroid_similarity_sha256\":\""+hash(similarity)+"\",\"class_energy_sha256\":\""+hash(energy)+"\"");
        System.out.print(",\"denominator_sha256\":\""+hash(denominator)+"\",\"updated_membership_sha256\":\""+hash(memberships)+"\",\"objective_bits\":"+Float.floatToIntBits(objective)+"}");
    }
}
