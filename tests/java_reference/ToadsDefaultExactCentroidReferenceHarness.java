import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected default 10obj exact computeCentroids composition after thinning. */
public final class ToadsDefaultExactCentroidReferenceHarness {
    private static final int N=36,C=11,BORDER=10;
    private static final float ZERO=1.0e-10f, SMOOTH=0.1f;
    private static final float[][] CENTER=new float[1][C];
    private static final float[] EXPONENT={0.80f,1.00f,1.00f};
    private static final String SEG_SHA="63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String ATLAS_SHA="cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939";
    private static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private final float[][][][] images=new float[N][N][N][1],field=new float[N][N][N][1],mems=new float[N][N][N][C];
    private final byte[][][] seg=new byte[N][N][N];
    private final float[] numerator=new float[C],denominator=new float[C];
    private final float[] updated=new float[C];

    private static int label(int x,int y,int z){return Math.floorMod(x+2*y+3*z+(y*z)%3+(x*y)%2,C);}
    private static float intensity(int x,int y,int z,int k){int delta=(x*13+y*7+z*3)%17-8;return (k+1)/12.0f+delta/1700.0f;}
    private static float fieldValue(int x,int y,int z){return 1.0f+((x+3*y+5*z)%9)/100.0f;}

    private void initialize(){
        for(int k=0;k<C;k++)CENTER[0][k]=(k+1)/12.0f;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int current=label(x,y,z);seg[x][y][z]=(byte)current;
            images[x][y][z][0]=intensity(x,y,z,current);field[x][y][z][0]=fieldValue(x,y,z);
        }
        // Selected exact initial membership method from MultipleToadPreprocess.
        for(int x=BORDER;x<N-BORDER-1;x++)for(int y=BORDER;y<N-BORDER-1;y++)for(int z=BORDER;z<N-BORDER-1;z++){
            float den=0.0f;
            for(int k=0;k<C;k++){
                float difference=images[x][y][z][0]-CENTER[0][k];float distance=difference*difference;
                if(distance>1.0e-20f)distance=1.0f/distance;else distance=1.0e20f;
                mems[x][y][z][k]=distance;den+=distance;
            }
            if(den>0.0f)for(int k=0;k<C;k++)mems[x][y][z][k]=mems[x][y][z][k]/den;
        }
    }
    private void update(){
        PowerTable[] table=new PowerTable[3];
        for(int category=0;category<3;category++)table[category]=new PowerTable(0.0f,4.0f,0.000001f,EXPONENT[category]);
        for(int k=0;k<C;k++){
            float num=0.0f,den=0.0f;
            int category=k==0?-1:k<=2?0:k<=7?1:2;
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
                float img=field[x][y][z][0]*images[x][y][z][0];
                float mem=mems[x][y][z][k];float weight;
                if(category<0)weight=mem*mem;
                else {
                    float residual=img-CENTER[0][k];float distance=residual*residual;float factor=1.0f;
                    if(distance>ZERO)factor=(float)(table[category].lookup((double)distance,(double)EXPONENT[category])/distance);
                    weight=factor*mem*mem;
                }
                num+=weight*img;den+=weight;
            }
            numerator[k]=num;denominator[k]=den;
            if(den>0.0f)updated[k]=(1.0f-SMOOTH)*num/den+SMOOTH*CENTER[0][k];
            else updated[k]=CENTER[0][k];
        }
    }
    private static String bits(float[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(values[i]));}return s.append(']').toString();}
    private static String values(float[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.toString(values[i]));}return s.append(']').toString();}
    public static void main(String[] args){
        ToadsDefaultExactCentroidReferenceHarness h=new ToadsDefaultExactCentroidReferenceHarness();h.initialize();h.update();
        System.out.println("{\"source\":\"VandyCRUISE MultipleToadSegmentation.computeExactCentroids default MASK+optimized 10obj dispatch\","
          +"\"segmentation_source_sha256\":\""+SEG_SHA+"\",\"preprocess_source_sha256\":\""+PREPROCESS_SHA+"\",\"atlas_descriptor_sha256\":\""+ATLAS_SHA+"\",\"power_table_source_sha256\":\""+POWER_SHA+"\","
          +"\"shape_xyz_channels\":[36,36,36,1],\"classes\":11,\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],"
          +"\"category_exponents\":[0.8,1.0,1.0],\"smoothness\":0.1,\"centroids\":"+values(CENTER[0])+",\"numerator_bits\":"+bits(h.numerator)+",\"denominator_bits\":"+bits(h.denominator)+",\"updated_centroids_bits\":"+bits(h.updated)+"}");
    }
}
