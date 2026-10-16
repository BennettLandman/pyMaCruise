/*
 * Selected AtlasEMPreprocess initialization methods, copied from the pinned
 * source class. This intentionally omits image registration/cropping.
 */
public final class AtlasEMPreprocessReferenceHarness {
    static final String SOURCE_SHA="acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b";
    static final int NX=4, NY=3, NZ=3, NC=2, K=4;
    static final float[][][][] images=new float[NX][NY][NZ][NC];
    static final float[][] centroid=new float[NC][K];
    static final float ZERO=1e-20f, INF=1e20f;
    static void init(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<NC;c++)
            images[x][y][z][c]=(float)(x*13+y*5+z*3+c*7+1)/17.0f;
        float[] low={0.25f,0.5f}, high={2.25f,4.5f};
        float[][] priors={{0.0f,0.3f,0.3f,0.8f},{0.1f,0.2f,0.4f,0.7f}};
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++)
            centroid[c][k]=low[c]/(high[c]-low[c])+priors[c][k];
    }
    static float[][][][] memberships(){
        float[][][][] out=new float[NX][NY][NZ][K];
        for(int x=0;x<NX-1;x++)for(int y=0;y<NY-1;y++)for(int z=0;z<NZ-1;z++){
            float den=0.0f, num;
            float[] grouping=new float[K];
            for(int k=0;k<K;k++){
                num=0.0f;
                for(int c=0;c<NC;c++)num+=(images[x][y][z][c]-centroid[c][k])*(images[x][y][z][c]-centroid[c][k]);
                if(num>ZERO)num=1.0f/num;else num=INF;
                out[x][y][z][k]=num; den+=num;
                grouping[k]=0.0f;
                for(int m=0;m<K;m++)if(centroid[0][m]==centroid[0][k])grouping[k]++;
            }
            if(den>0.0f)for(int k=0;k<K;k++)out[x][y][z][k]=grouping[k]*out[x][y][z][k]/den;
        }
        return out;
    }
    static String vector(float[] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=Float.floatToRawIntBits(a[i]);}return s+"]";}
    static String matrix(float[][] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=vector(a[i]);}return s+"]";}
    static String tensor(float[][][][] a){String s="[";boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){if(!first)s+=",";first=false;s+=Float.floatToRawIntBits(a[x][y][z][k]);}return s+"]";}
    public static void main(String[] args){init();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"centroid\":"+matrix(centroid)+",\"memberships\":"+tensor(memberships())+"}");}
}
