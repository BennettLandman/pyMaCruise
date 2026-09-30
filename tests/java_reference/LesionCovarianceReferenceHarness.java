/** Selected-method excerpt of VandyCRUISE LesionToadSegmentation.computCovariance. */
public final class LesionCovarianceReferenceHarness {
    private static final int NX=5, NY=4, NZ=3, NC=3, C=6;
    private static final short[][][] seg = new short[NX][NY][NZ];
    private static final float[][][][] mem = new float[NX][NY][NZ][C];
    private static final float[][][][] image = new float[NX][NY][NZ][NC];
    private static final float[][][][] field = new float[NX][NY][NZ][NC];
    private static final float[][] centroid = new float[NC][C];
    private static final float[] lesionCentroid = {0.25f, 0.5625f, 0.8125f};
    private static final float[][][] lesion = new float[NX][NY][NZ];
    private static final String[] type = {"object", "object", "lesion", "blackhole", "outlier", "object"};

    private static int bits(float value) { return Float.floatToRawIntBits(value); }

    private static final class Result {
        float[][][] classes = new float[C][NC][NC];
        float[][] lesion = new float[NC][NC];
    }

    // This follows the selected method's class -> x/y/z -> channel-pair loops,
    // including its shared lesion accumulator and normalization at LESION class.
    private static Result compute(boolean useField) {
        Result out = new Result();
        for (int c1=0; c1<NC; c1++) for (int c2=0; c2<NC; c2++) out.lesion[c1][c2]=0.0f;
        for (int k=0; k<C; k++) {
            for (int c1=0; c1<NC; c1++) for (int c2=0; c2<NC; c2++) out.classes[k][c1][c2]=0.0f;
            float den=0.0f, denL=0.0f;
            for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) if (seg[x][y][z] == k) {
                boolean toLesion = type[k].equals("lesion") && mem[x][y][z][k] < lesion[x][y][z];
                float m = toLesion ? lesion[x][y][z] : mem[x][y][z][k];
                float weight = m*m + 1e-30f;
                for (int c1=0; c1<NC; c1++) for (int c2=0; c2<NC; c2++) {
                    float sample1 = useField ? field[x][y][z][c1]*image[x][y][z][c1] : image[x][y][z][c1];
                    float sample2 = useField ? field[x][y][z][c2]*image[x][y][z][c2] : image[x][y][z][c2];
                    float mean1 = toLesion ? lesionCentroid[c1] : centroid[c1][k];
                    float mean2 = toLesion ? lesionCentroid[c2] : centroid[c2][k];
                    if (toLesion) out.lesion[c1][c2] += weight*(sample1-mean1)*(sample2-mean2);
                    else out.classes[k][c1][c2] += weight*(sample1-mean1)*(sample2-mean2);
                }
                if (toLesion) denL += weight;
                else den += weight;
            }
            for (int c1=0; c1<NC; c1++) for (int c2=0; c2<NC; c2++) {
                if (den != 0.0f) out.classes[k][c1][c2] /= den;
                else out.classes[k][c1][c2] = 0.0f;
                if (type[k].equals("lesion")) {
                    if (denL != 0.0f) out.lesion[c1][c2] /= denL;
                    else out.lesion[c1][c2] = 0.0f;
                }
            }
        }
        return out;
    }

    private static void init() {
        for (int c=0; c<NC; c++) for (int k=0; k<C; k++) centroid[c][k] = ((c*3+k*2+1)%13)*0.0625f;
        for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) {
            seg[x][y][z] = (short)((x+2*y+z)%5); // class 5 intentionally absent
            for (int k=0; k<C; k++) {
                int code=(x*3+y*5+z*7+k*2)%9;
                mem[x][y][z][k]=(code+1)*0.0625f;
                if (k==seg[x][y][z] && (x+2*y+z)%11==0) mem[x][y][z][k]=0.0f;
            }
            for (int c=0; c<NC; c++) {
                image[x][y][z][c]=((x*11+y*7+z*3+c*5)%31)*0.0625f;
                field[x][y][z][c]=1.0f+((x*2+y*3+z*4+c)%5)*0.0625f;
            }
            float l=((x+2*y+z)%7+1)*0.0625f;
            if (seg[x][y][z]==2 && (x+y+z)%3==0) l=mem[x][y][z][2]; // strict-comparison tie
            lesion[x][y][z]=l;
        }
    }

    private static String matrix(float[][] a) {
        StringBuilder b=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('[');for(int j=0;j<a[i].length;j++){if(j>0)b.append(',');b.append(bits(a[i][j]));}b.append(']');}
        return b.append(']').toString();
    }
    private static String cube(float[][][] a) {
        StringBuilder b=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(matrix(a[i]));}
        return b.append(']').toString();
    }
    private static String json(Result r) { return "{\"class_bits\":"+cube(r.classes)+",\"lesion_bits\":"+matrix(r.lesion)+"}"; }
    public static void main(String[] args) {
        init();
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"no_field\":"+json(compute(false))+",\"with_field\":"+json(compute(true))+"}");
    }
}
