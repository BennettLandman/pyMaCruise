import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships RAYLEIGH data-term excerpt. */
public final class ToadsRayleighIntensityEnergyReferenceHarness {
    private static final String SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=11, NY=9, NZ=7, CHANNELS=3, CLASSES=4;
    private static final float[][] CENTROIDS={
        {0.25f,0.40f,0.0f,-0.40f},
        {0.60f,0.30f,0.70f,0.50f},
        {0.80f,0.55f,0.45f,0.65f}
    };
    private static final float[] MODALITY={0.73f,0.41f,0.89f};
    private static final float[][] CLASS_VARIANCE={
        {0.61f,1.13f,0.77f,0.52f},
        {1.21f,0.68f,0.93f,1.07f},
        {0.84f,1.31f,0.59f,0.72f}
    };
    private static final float[] IMAGE_VARIANCE={0.57f,0.32f,0.76f};
    private static float[][][][] images, field, bias;

    private static byte[] readGzip(String path) throws Exception {
        ByteArrayOutputStream out=new ByteArrayOutputStream();
        try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){
            byte[] buffer=new byte[8192]; int n;
            while((n=in.read(buffer))>=0) out.write(buffer,0,n);
        }
        return out.toByteArray();
    }
    private static float[][][][] read4d(String path) throws Exception {
        byte[] raw=readGzip(path); ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348 || b.getShort(70)!=16 || b.getShort(72)!=32 || b.getShort(40)!=4 ||
                b.getShort(42)!=NX || b.getShort(44)!=NY || b.getShort(46)!=NZ || b.getShort(48)!=CHANNELS)
            throw new IllegalArgumentException("unexpected float32 RAYLEIGH fixture header");
        int offset=(int)b.getFloat(108); float[][][][] data=new float[NX][NY][NZ][CHANNELS];
        for(int c=0;c<CHANNELS;c++) for(int z=0;z<NZ;z++) for(int y=0;y<NY;y++) for(int x=0;x<NX;x++){
            int i=(((c*NZ+z)*NY+y)*NX+x); data[x][y][z][c]=b.getFloat(offset+4*i);
        }
        return data;
    }
    private static float rayleighDistance(float val,float sig){
        final float logpi2mpi2=(float)Math.log(Math.PI/2.0f)-(float)Math.PI/2.0f;
        if(val<=0) return 1.0f;
        else return val*val/sig*sig-2.0f*(float)Math.log(val/sig)+logpi2mpi2;
    }
    private static int[] energyBits(boolean useField,boolean useBias){
        float[][][][] energy=new float[NX][NY][NZ][CLASSES];
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++){
            for(int k=0;k<CLASSES;k++){
                float num=0.0f;
                for(int c=0;c<CHANNELS;c++){
                    float value=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];
                    if(useBias) value-=bias[x][y][z][c];
                    num+=MODALITY[c]*CLASS_VARIANCE[c][k]*IMAGE_VARIANCE[c]*rayleighDistance(value,CENTROIDS[c][k]);
                }
                energy[x][y][z][k]=num;
            }
        }
        int[] bits=new int[NX*NY*NZ*CLASSES]; int p=0;
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int k=0;k<CLASSES;k++)
            bits[p++]=Float.floatToIntBits(energy[x][y][z][k]);
        return bits;
    }
    private static void printInts(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    private static void printFloatBits(StringBuilder s,float[] a){int[] bits=new int[a.length];for(int i=0;i<a.length;i++)bits[i]=Float.floatToIntBits(a[i]);printInts(s,bits);}
    private static void printMatrixBits(StringBuilder s,float[][] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');printFloatBits(s,a[i]);}s.append(']');}

    public static void main(String[] args)throws Exception{
        images=read4d(args[0]);field=read4d(args[1]);bias=read4d(args[2]);
        StringBuilder out=new StringBuilder();
        out.append("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[11,9,7,4],\"centroid_bits\":");
        printMatrixBits(out,CENTROIDS);out.append(",\"modality_bits\":");printFloatBits(out,MODALITY);
        out.append(",\"class_variance_bits\":");printMatrixBits(out,CLASS_VARIANCE);
        out.append(",\"image_variance_bits\":");printFloatBits(out,IMAGE_VARIANCE);
        out.append(",\"cases\":[");
        boolean[][] modes={{false,false},{false,true},{true,false},{true,true}};
        for(int i=0;i<modes.length;i++){
            if(i>0)out.append(',');
            out.append("{\"field\":").append(modes[i][0]).append(",\"bias\":").append(modes[i][1]).append(",\"energy_bits\":");
            printInts(out,energyBits(modes[i][0],modes[i][1]));out.append('}');
        }
        out.append("]}");System.out.println(out);
    }
}
