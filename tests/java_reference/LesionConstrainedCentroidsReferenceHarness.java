/** Selected-method excerpt of LesionToadSegmentation.computeConstrainedCentroids. */
public final class LesionConstrainedCentroidsReferenceHarness {
    private static final int NX=5, NY=4, NZ=3, NC=2, C=4;
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC];
    private static final float[][][][] field=new float[NX][NY][NZ][NC];
    private static final float[][] INITIAL={{0.15f,0.55f,1.05f,1.65f},{0.30f,0.70f,1.15f,1.80f}};

    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[][] update(boolean useField){
        float[][] centroid=new float[NC][C];for(int c=0;c<NC;c++)centroid[c]=INITIAL[c].clone();
        float factor=0.45f;
        for(int c=0;c<NC;c++){
            float prev=-1.0e10f;
            for(int k=0;k<C;k++){
                float num=0.0f,den=0.0f;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                    float m=mem[x][y][z][k];float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                    num+=m*m*sample;den+=m*m;
                }
                float next=den>0.0f?num/den:0.0f;
                if(next<centroid[c][k]-factor*(centroid[c][k]-prev))next=centroid[c][k]-factor*(centroid[c][k]-prev);
                else if(k+1<C&&next>centroid[c][k]+factor*(centroid[c][k+1]-centroid[c][k]))next=centroid[c][k]+factor*(centroid[c][k+1]-centroid[c][k]);
                prev=centroid[c][k];centroid[c][k]=next;
            }
        }
        return centroid;
    }
    private static String bitsJson(float[][] values){StringBuilder b=new StringBuilder("[");for(int c=0;c<NC;c++){if(c>0)b.append(',');b.append('[');for(int k=0;k<C;k++){if(k>0)b.append(',');b.append(bits(values[c][k]));}b.append(']');}return b.append(']').toString();}
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int label=(x+2*y+z)%3;seg[x][y][z]=(short)label;
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
        }
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"without_field_bits\":"+bitsJson(update(false))+",\"with_field_bits\":"+bitsJson(update(true))+"}");
    }
}
