/**
 * Independent excerpt harness for VandyCRUISE ordinary OBJECT exact centroids.
 * It records the x/y/z scan, optional field correction, boundary PSF-ridge
 * arithmetic, and active AlgorithmMultipleToads damped mode (smoothness .1).
 * This is not a port of MultipleToadSegmentation's dispatcher.
 */
public final class ToadsObjectCentroidsReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String WRAPPER_SHA =
        "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final int NX=5, NY=4, NZ=4, CHANNELS=2, CLASSES=3;
    private static final float SMOOTH=0.1f, RIDGE=0.23f;
    private static float[][][][] image = new float[NX][NY][NZ][CHANNELS];
    private static float[][][][] field = new float[NX][NY][NZ][CHANNELS];
    private static float[][][][] mem = new float[NX][NY][NZ][CLASSES];
    private static short[][][] seg = new short[NX][NY][NZ];
    private static float[][] initial = {{0.31f, 0.72f, 1.19f}, {0.43f, 0.88f, 1.37f}};
    private static int bits(float value) { return Float.floatToRawIntBits(value); }

    private static void setup() {
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            seg[x][y][z]=(short)((x+2*y+z)%3);
            for(int c=0;c<CHANNELS;c++) {
                image[x][y][z][c]=0.17f+((x*7+y*5+z*3+c*11)%19)*0.071f;
                field[x][y][z][c]=0.73f+((x*3+y*2+z*5+c)%7)*0.047f;
            }
            for(int q=0;q<CLASSES;q++) {
                mem[x][y][z][q]=q==2?0.0f:0.09f+((x*2+y*3+z*5+q*7)%11)*0.063f;
            }
        }
    }

    private static boolean isBoundary(int x,int y,int z) {
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
            if(seg[x+i][y+j][z+l]!=seg[x][y][z]) return true;
        return false;
    }

    private static float[][] update(boolean useField, boolean useRidge) {
        float[][] result=new float[CHANNELS][CLASSES];
        for(int c=0;c<CHANNELS;c++)result[c]=initial[c].clone();
        for(int c=0;c<CHANNELS;c++) {
            float[] num=new float[CLASSES],den=new float[CLASSES];
            for(int q=0;q<CLASSES;q++)for(short x=1;x<NX-1;x++)for(short y=1;y<NY-1;y++)for(short z=1;z<NZ-1;z++) {
                float m=mem[x][y][z][q];
                float img;
                if(useField)img=field[x][y][z][c]*image[x][y][z][c];
                else img=image[x][y][z][c];
                if(useRidge && isBoundary(x,y,z)) {
                    img+=RIDGE*(28.0f/27.0f)*image[x][y][z][c];
                    for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                        img-=RIDGE*(1.0f/27.0f)*image[x+i][y+j][z+l][c];
                }
                num[q]+=m*m*1.0f*img;
                den[q]+=m*m*1.0f*1.0f;
            }
            for(int q=0;q<CLASSES;q++)if(den[q]>0.0f)
                result[c][q]=(1.0f-SMOOTH)*num[q]/den[q]+SMOOTH*result[c][q];
        }
        return result;
    }

    private static void append(StringBuilder out,float[][] a) {
        out.append('[');
        for(int c=0;c<CHANNELS;c++) {
            if(c>0)out.append(',');out.append('[');
            for(int q=0;q<CLASSES;q++){if(q>0)out.append(',');out.append(bits(a[c][q]));}
            out.append(']');
        }
        out.append(']');
    }
    public static void main(String[] args) {
        setup();StringBuilder b=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"wrapper_sha256\":\"").append(WRAPPER_SHA).append("\"");
        b.append(",\"raw\":");append(b,update(false,false));
        b.append(",\"field\":");append(b,update(true,false));
        b.append(",\"ridge\":");append(b,update(false,true));
        b.append(",\"field_ridge\":");append(b,update(true,true));
        System.out.println(b.append('}'));
    }
}
