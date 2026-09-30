import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships OBJECT/MASK/OUTMASK source composition. */
public final class ToadsOutmaskMembershipPassReferenceHarness {
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=9,NY=8,NZ=7,CHANNELS=2,CLASSES=3,OUTMASK=1;
    private static final String[] TYPES={"MASK","OUTMASK","OBJECT"};
    private static final float[][] CENTROIDS={{0.23f,0.31f,0.42f},{0.38f,0.27f,0.51f}};
    private static final float[] MODALITY={0.73f,0.46f};
    private static final float[][] VARIANCE={{0.81f,1.17f,0.69f},{1.24f,0.77f,0.95f}};
    private static final float[] IMAGE_VARIANCE={0.61f,0.34f};
    private static final float[] SMOOTH_PRIOR={0.65f,0.91f,1.13f};
    private static final float[] ATLAS_PRIOR={0.20f,0.33f,0.47f};
    private static final float[][] SIMILARITY={{0,0.17f,0.24f},{0.13f,0,0.19f},{0.22f,0.16f,0}};
    private static final float OUTLIER_VALUE=1.37f,SMOOTHING=0.043f,PRIOR_FACTOR=0.64f,PRIOR_COEFFICIENT=0.16f;
    private static final float ZERO=1e-10f,INF=1e10f;
    private static float[][][][] images,field,bias,initial,shapes;
    private static float[][][] initialOutlier;
    private static short[][][][] labels;

    private static byte[] unzip(String path)throws Exception{
        ByteArrayOutputStream out=new ByteArrayOutputStream();try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] a=new byte[8192];int n;while((n=in.read(a))>=0)out.write(a,0,n);}return out.toByteArray();
    }
    private static float[][][][] readFloat4(String path,int channels)throws Exception{
        byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=channels)throw new IllegalArgumentException("unexpected float NIfTI fixture");
        int offset=(int)b.getFloat(108);float[][][][] a=new float[NX][NY][NZ][channels];
        for(int c=0;c<channels;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getFloat(offset+4*(((c*NZ+z)*NY+y)*NX+x));return a;
    }
    private static float[][][] readFloat3(String path)throws Exception{
        byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=3||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ)throw new IllegalArgumentException("unexpected float3 NIfTI fixture");
        int offset=(int)b.getFloat(108);float[][][] a=new float[NX][NY][NZ];for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z]=b.getFloat(offset+4*((z*NY+y)*NX+x));return a;
    }
    private static short[][][][] readShort4(String path,int channels)throws Exception{
        byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=channels)throw new IllegalArgumentException("unexpected int16 NIfTI fixture");
        int offset=(int)b.getFloat(108);short[][][][] a=new short[NX][NY][NZ][channels];for(int c=0;c<channels;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getShort(offset+2*(((c*NZ+z)*NY+y)*NX+x));return a;
    }
    private static float maskDistance(float image){return image==0.0f?0.0f:1.0f;}
    private static float[][][][] dataEnergy(boolean useField,boolean useBias,float[][][] outlierEnergy){
        float[][][][] energy=new float[NX][NY][NZ][CLASSES];
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<CLASSES;k++){
            float num=0.0f,numO=0.0f;
            if(TYPES[k].equals("OBJECT"))for(int c=0;c<CHANNELS;c++){
                float value=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];
                if(useBias)value-=bias[x][y][z][c];float residual=value-CENTROIDS[c][k];float weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num+=weight*residual*residual;
            }
            else if(TYPES[k].equals("MASK")||TYPES[k].equals("OUTMASK"))for(int c=0;c<CHANNELS;c++){
                float weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num+=weight*maskDistance(images[x][y][z][c]);
                if(TYPES[k].equals("OUTMASK")){numO+=weight*OUTLIER_VALUE;outlierEnergy[x][y][z]=numO;}
            }
            energy[x][y][z][k]=num;
        }
        return energy;
    }
    private static float shapeFactor(float s){return Math.min(1.0f,2.0f*PRIOR_FACTOR*s);}
    private static float neighborFactor(float s){return Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*s,1.0f));}
    private static float[][][][] update(boolean useField,boolean useBias,float[][][] outlierEnergy,float[] total){
        float[][][][] energy=dataEnergy(useField,useBias,outlierEnergy),mem=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)mem[x][y][z][k]=initial[x][y][z][k];
        float[][][] out=new float[NX][NY][NZ];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)out[x][y][z]=initialOutlier[x][y][z];
        float objective=0.0f;int[][] dirs={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0.0f;float[] objAt=new float[CLASSES];
            for(int k=0;k<CLASSES;k++){
                float num=energy[x][y][z][k],numO=outlierEnergy[x][y][z],obj=num,currentOut=out[x][y][z];
                if(k==OUTMASK)obj=(1.0f-currentOut)*num+currentOut*numO;
                if(SMOOTHING>0.0f){float ngb=0.0f;for(int[] d:dirs)for(int m=0;m<CLASSES;m++)if(m!=k){float v=mem[x+d[0]][y+d[1]][z+d[2]][m];ngb+=v*v;}
                    float increment=SMOOTH_PRIOR[k]*SMOOTHING*ngb;num+=increment;obj+=0.5f*SMOOTH_PRIOR[k]*SMOOTHING*ngb;
                    float numL=0.0f;
                    // Preserve the selected source's unbraced nested-if/else
                    // chain: OUTMASK never enters the auxiliary smoothing arm.
                    if(TYPES[k].equals("LESION"))if(false)numL+=increment;else if(TYPES[k].equals("OUTLIER"))numO+=increment;else if(TYPES[k].equals("OUTMASK"))numO+=increment;
                }
                if(PRIOR_COEFFICIENT>0.0f){float inverse=0.0f;for(int m=0;m<CLASSES;m++)if(m!=k)inverse+=ATLAS_PRIOR[m]*SIMILARITY[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];
                    num+=PRIOR_COEFFICIENT*inverse;if(k==OUTMASK)numO+=PRIOR_COEFFICIENT*inverse;obj+=0.5f*PRIOR_COEFFICIENT*inverse;}
                float reciprocal=num>ZERO?1.0f/num:INF;boolean direct=labels[x][y][z][0]==k||labels[x][y][z][1]==k;boolean near=labels[x][y][z][2]==k||labels[x][y][z][3]==k;
                float factor=direct?1.0f:(near?neighborFactor(shapes[x][y][z][k]):shapeFactor(shapes[x][y][z][k]));
                float value;
                if(k==OUTMASK){float reciprocalO=numO>ZERO?1.0f/numO:INF;value=factor*((1.0f-currentOut)*reciprocal+currentOut*reciprocalO);out[x][y][z]=(factor*currentOut)*reciprocalO;}
                else value=factor*reciprocal;
                mem[x][y][z][k]=value;den+=value;objAt[k]=obj;
            }
            if(den>0.0f)for(int k=0;k<CLASSES;k++){mem[x][y][z][k]=mem[x][y][z][k]/den;if(k==OUTMASK)out[x][y][z]=out[x][y][z]/den;objective+=objAt[k]*mem[x][y][z][k]*mem[x][y][z][k];}
        }
        total[0]=objective;lastOutlier=out;return mem;
    }
    private static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    private static int[] flatten(float[][][][] a){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(a[x][y][z][k]);return b;}
    private static int[] flatten(float[][][] a){int[] b=new int[NX*NY*NZ];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)b[p++]=Float.floatToIntBits(a[x][y][z]);return b;}
    private static int[] bits(float[] a){int[] b=new int[a.length];for(int i=0;i<a.length;i++)b[i]=Float.floatToIntBits(a[i]);return b;}
    private static int[] bits(float[][] a){int[] b=new int[a.length*a[0].length];int p=0;for(float[] row:a)for(float v:row)b[p++]=Float.floatToIntBits(v);return b;}

    public static void main(String[] args)throws Exception{
        images=readFloat4(args[0],CHANNELS);field=readFloat4(args[1],CHANNELS);bias=readFloat4(args[2],CHANNELS);initial=readFloat4(args[3],CLASSES);shapes=readFloat4(args[4],CLASSES);labels=readShort4(args[5],4);initialOutlier=readFloat3(args[6]);
        StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[9,8,7,3],\"class_types\":[\"MASK\",\"OUTMASK\",\"OBJECT\"]");
        s.append(",\"centroid_bits\":");append(s,bits(CENTROIDS));s.append(",\"modality_bits\":");append(s,bits(MODALITY));s.append(",\"variance_bits\":");append(s,bits(VARIANCE));s.append(",\"image_variance_bits\":");append(s,bits(IMAGE_VARIANCE));s.append(",\"smoothing_prior_bits\":");append(s,bits(SMOOTH_PRIOR));s.append(",\"atlas_prior_bits\":");append(s,bits(ATLAS_PRIOR));
        int[] sim=new int[CLASSES*CLASSES];int p=0;for(float[] row:SIMILARITY)for(float v:row)sim[p++]=Float.floatToIntBits(v);s.append(",\"similarity_bits\":");append(s,sim);
        s.append(",\"outlier_value_bits\":").append(Float.floatToIntBits(OUTLIER_VALUE)).append(",\"smoothing\":").append(SMOOTHING).append(",\"prior_factor\":").append(PRIOR_FACTOR).append(",\"prior_coefficient\":").append(PRIOR_COEFFICIENT).append(",\"cases\":[");
        boolean[][] modes={{false,false},{false,true},{true,false},{true,true}};
        for(int i=0;i<modes.length;i++){if(i>0)s.append(',');float[][][] outEnergy=new float[NX][NY][NZ];float[] objective={0};float[][][][] mem=update(modes[i][0],modes[i][1],outEnergy,objective);
            float[][][] classOut=new float[NX][NY][NZ];for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)classOut[x][y][z]=outEnergy[x][y][z];
            s.append("{\"field\":").append(modes[i][0]).append(",\"bias\":").append(modes[i][1]).append(",\"class_energy_bits\":");append(s,flatten(dataEnergy(modes[i][0],modes[i][1],new float[NX][NY][NZ])));
            s.append(",\"outlier_energy_bits\":");append(s,flatten(classOut));s.append(",\"membership_bits\":");append(s,flatten(mem));s.append(",\"outlier_membership_bits\":");append(s,flatten(lastOutlier));s.append(",\"objective_bits\":").append(Float.floatToIntBits(objective[0])).append('}');}
        s.append("]}");System.out.println(s);
    }
    private static float[][][] lastOutlier;
}
