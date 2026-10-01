/** Selected full-grid computeExactVariances method excerpt for default 10obj. */
public final class ToadsDefaultExactVarianceReferenceHarness {
    private static final int N=36,C=11,BORDER=10;
    private static final float[][] CENTER=new float[1][C];
    private static final int[] REFRESHED_BITS={1035019062,1043302233,1048826665,1051777520,1054717732,1057274465,1058720888,1060177692,1061602345,1063073864,1064536744};
    private final float[][][][] images=new float[N][N][N][1],field=new float[N][N][N][1],mems=new float[N][N][N][C];
    private final byte[][][] seg=new byte[N][N][N];
    private final float[] numerator=new float[C],denominator=new float[C],deviation=new float[C];

    private static int label(int x,int y,int z){return Math.floorMod(x+2*y+3*z+(y*z)%3+(x*y)%2,C);}
    private static float intensity(int x,int y,int z,int k){int delta=(x*13+y*7+z*3)%17-8;return (k+1)/12.0f+delta/1700.0f;}
    private static float fieldValue(int x,int y,int z){return 1.0f+((x+3*y+5*z)%9)/100.0f;}
    private void initialize(){
        for(int k=0;k<C;k++)CENTER[0][k]=(k+1)/12.0f;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int current=label(x,y,z);seg[x][y][z]=(byte)current;images[x][y][z][0]=intensity(x,y,z,current);field[x][y][z][0]=fieldValue(x,y,z);
        }
        // Source initialExactMemberships on the selected ten-voxel crop.
        for(int x=BORDER;x<N-BORDER-1;x++)for(int y=BORDER;y<N-BORDER-1;y++)for(int z=BORDER;z<N-BORDER-1;z++){
            float den=0.0f;
            for(int k=0;k<C;k++){
                float difference=images[x][y][z][0]-CENTER[0][k];float distance=difference*difference;
                if(distance>1.0e-20f)distance=1.0f/distance;else distance=1.0e20f;
                mems[x][y][z][k]=distance;den+=distance;
            }
            if(den>0.0f)for(int k=0;k<C;k++)mems[x][y][z][k]=mems[x][y][z][k]/den;
        }
        for(int k=0;k<C;k++)CENTER[0][k]=Float.intBitsToFloat(REFRESHED_BITS[k]);
    }
    private void update(){
        // MultipleToadSegmentation.computeExactVariances, non-LESION branch.
        for(int k=0;k<C;k++){
            float num=0.0f,den=0.0f;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)if(seg[x][y][z]==k){
                float mem=mems[x][y][z][k];float squared=mem*mem;float weight=squared+1.0e-30f;
                float corrected=field[x][y][z][0]*images[x][y][z][0];float diff=corrected-CENTER[0][k];
                num+=weight*diff*diff;den+=weight;
            }
            numerator[k]=num;denominator[k]=den;
            deviation[k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
    }
    private static String bits(float[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(values[i]));}return s.append(']').toString();}
    public static void main(String[] args){
        ToadsDefaultExactVarianceReferenceHarness h=new ToadsDefaultExactVarianceReferenceHarness();h.initialize();h.update();
        System.out.println("{\"source\":\"VandyCRUISE MultipleToadSegmentation.computeExactVariances default MASK+optimized 10obj dispatch\","
          +"\"segmentation_source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\","
          +"\"preprocess_source_sha256\":\"31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7\","
          +"\"atlas_descriptor_sha256\":\"cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939\","
          +"\"shape_xyz_channels\":[36,36,36,1],\"classes\":11,\"centroids_bits\":"+bits(CENTER[0])+",\"numerator_bits\":"+bits(h.numerator)+",\"denominator_bits\":"+bits(h.denominator)+",\"standard_deviation_bits\":"+bits(h.deviation)+"}");
    }
}
