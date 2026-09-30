import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships RAYLEIGH/data/prior/update composition. */
public final class ToadsRayleighMembershipCompositionReferenceHarness {
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=9,NY=8,NZ=7,CHANNELS=3,CLASSES=3;
    private static final float[][] CENTROIDS={{0.27f,0.34f,0.42f},{0.31f,0.39f,0.29f},{0.46f,0.36f,0.51f}};
    private static final float[] MODALITY={0.71f,0.53f,0.36f};
    private static final float[][] VARIANCE={{0.82f,1.13f,0.67f},{1.21f,0.74f,0.93f},{0.58f,1.32f,0.77f}};
    private static final float[] IMAGE_VARIANCE={0.63f,0.41f,0.28f};
    private static final float[] SMOOTH_PRIOR={0.6f,0.9f,1.1f};
    private static final float[] ATLAS_PRIOR={0.20f,0.35f,0.45f};
    private static final float[][] SIMILARITY={{0,0.17f,0.21f},{0.13f,0,0.19f},{0.23f,0.16f,0}};
    private static final float SMOOTHING=0.032f,PRIOR_FACTOR=0.68f,PRIOR_COEFFICIENT=0.14f;
    private static final float ZERO=1e-10f,INF=1e10f;
    private static float[][][][] images,field,bias,initial,shapes;
    private static short[][][][] labels;
    private static final float LOGPI2MPI2=(float)Math.log(Math.PI/2.0f)-(float)Math.PI/2.0f;

    private static byte[] gunzip(String path)throws Exception{
        ByteArrayOutputStream out=new ByteArrayOutputStream();
        try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] buf=new byte[8192];int n;while((n=in.read(buf))>=0)out.write(buf,0,n);}
        return out.toByteArray();
    }
    private static float[][][][] readFloat4(String path,int channels)throws Exception{
        byte[] raw=gunzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=channels)
            throw new IllegalArgumentException("unexpected float NIfTI fixture");
        int offset=(int)b.getFloat(108);float[][][][] a=new float[NX][NY][NZ][channels];
        for(int c=0;c<channels;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getFloat(offset+4*(((c*NZ+z)*NY+y)*NX+x));
        return a;
    }
    private static short[][][][] readShort4(String path,int channels)throws Exception{
        byte[] raw=gunzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=channels)
            throw new IllegalArgumentException("unexpected int16 labels NIfTI fixture");
        int offset=(int)b.getFloat(108);short[][][][] a=new short[NX][NY][NZ][channels];
        for(int c=0;c<channels;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getShort(offset+2*(((c*NZ+z)*NY+y)*NX+x));
        return a;
    }
    private static float rayleigh(float val,float sig){if(val<=0)return 1.0f;return val*val/sig*sig-2.0f*(float)Math.log(val/sig)+LOGPI2MPI2;}
    private static float[][][][] energy(boolean useField,boolean useBias){
        float[][][][] e=new float[NX][NY][NZ][CLASSES];
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<CLASSES;k++){
            float num=0.0f;
            for(int c=0;c<CHANNELS;c++){
                float value=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];
                if(useBias)value-=bias[x][y][z][c];
                num+=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c]*rayleigh(value,CENTROIDS[c][k]);
            }
            e[x][y][z][k]=num;
        }
        return e;
    }
    private static float shapeFactor(float shape){return Math.min(1.0f,2.0f*PRIOR_FACTOR*shape);}
    private static float neighborFactor(float shape){return Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*shape,1.0f));}
    private static int[] run(boolean useField,boolean useBias){
        float[][][][] energy=energy(useField,useBias),mem=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)mem[x][y][z][k]=initial[x][y][z][k];
        float objective=0.0f;int[][] dirs={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0.0f;float[] objAt=new float[CLASSES];
            for(int k=0;k<CLASSES;k++){
                float num=energy[x][y][z][k],obj=num;
                float ngb=0.0f;
                if(SMOOTHING>0.0f){
                    for(int[] d:dirs)for(int m=0;m<CLASSES;m++)if(m!=k){float v=mem[x+d[0]][y+d[1]][z+d[2]][m];ngb+=v*v;}
                    num+=SMOOTH_PRIOR[k]*SMOOTHING*ngb;
                    obj+=0.5f*SMOOTH_PRIOR[k]*SMOOTHING*ngb;
                }
                float invprior=0.0f;
                for(int m=0;m<CLASSES;m++)if(m!=k){invprior+=ATLAS_PRIOR[m]*SIMILARITY[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];}
                if(PRIOR_COEFFICIENT>0.0f){num+=PRIOR_COEFFICIENT*invprior;obj+=0.5f*PRIOR_COEFFICIENT*invprior;}
                boolean direct=labels[x][y][z][0]==k||labels[x][y][z][1]==k;
                boolean near=labels[x][y][z][2]==k||labels[x][y][z][3]==k;
                float factor=direct?1.0f:(near?neighborFactor(shapes[x][y][z][k]):shapeFactor(shapes[x][y][z][k]));
                float value=num>ZERO?1.0f/num:INF;
                mem[x][y][z][k]=factor*value;den+=mem[x][y][z][k];objAt[k]=obj;
            }
            if(den>0.0f)for(int k=0;k<CLASSES;k++){
                mem[x][y][z][k]=mem[x][y][z][k]/den;
                float contribution=objAt[k]*mem[x][y][z][k];contribution*=mem[x][y][z][k];objective+=contribution;
            }
        }
        int n=NX*NY*NZ*CLASSES;int[] bits=new int[n+1];int p=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)bits[p++]=Float.floatToIntBits(mem[x][y][z][k]);
        bits[n]=Float.floatToIntBits(objective);return bits;
    }
    private static void printBits(StringBuilder s,int[] bits){s.append('[');for(int i=0;i<bits.length;i++){if(i>0)s.append(',');s.append(bits[i]);}s.append(']');}
    private static void printFloats(StringBuilder s,float[] a){int[] b=new int[a.length];for(int i=0;i<a.length;i++)b[i]=Float.floatToIntBits(a[i]);printBits(s,b);}
    private static void printMatrix(StringBuilder s,float[][] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');printFloats(s,a[i]);}s.append(']');}
    private static int[] flatten(float[][][][] a){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(a[x][y][z][k]);return b;}

    public static void main(String[] args)throws Exception{
        images=readFloat4(args[0],CHANNELS);field=readFloat4(args[1],CHANNELS);bias=readFloat4(args[2],CHANNELS);
        initial=readFloat4(args[3],CLASSES);shapes=readFloat4(args[4],CLASSES);labels=readShort4(args[5],4);
        boolean[][] modes={{false,false},{false,true},{true,false},{true,true}};
        StringBuilder out=new StringBuilder();out.append("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[9,8,7,3],\"centroid_bits\":");printMatrix(out,CENTROIDS);
        out.append(",\"modality_bits\":");printFloats(out,MODALITY);out.append(",\"variance_bits\":");printMatrix(out,VARIANCE);out.append(",\"image_variance_bits\":");printFloats(out,IMAGE_VARIANCE);
        out.append(",\"smoothing_prior_bits\":");printFloats(out,SMOOTH_PRIOR);out.append(",\"atlas_prior_bits\":");printFloats(out,ATLAS_PRIOR);out.append(",\"similarity_bits\":");printMatrix(out,SIMILARITY);
        out.append(",\"smoothing\":").append(SMOOTHING).append(",\"prior_factor\":").append(PRIOR_FACTOR).append(",\"prior_coefficient\":").append(PRIOR_COEFFICIENT).append(",\"cases\":[");
        for(int i=0;i<modes.length;i++){if(i>0)out.append(',');float[][][][] e=energy(modes[i][0],modes[i][1]);int[] result=run(modes[i][0],modes[i][1]);int[] membershipBits=java.util.Arrays.copyOf(result,result.length-1);
            out.append("{\"field\":").append(modes[i][0]).append(",\"bias\":").append(modes[i][1]).append(",\"energy_bits\":");printBits(out,flatten(e));
            out.append(",\"membership_bits\":");printBits(out,membershipBits);out.append(",\"objective_bits\":").append(result[result.length-1]).append('}');}
        out.append("]}");System.out.println(out);
    }
}
