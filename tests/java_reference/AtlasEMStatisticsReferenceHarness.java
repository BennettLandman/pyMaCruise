/*
 * Fixed-input AtlasEM statistics excerpt from AtlasEMSegmentation.java.
 * The loops and float32 accumulation order mirror computePriors,
 * computeCentriods, and computeCovariance. This avoids MIPAV dependencies and
 * does not model atlas registration or EM membership evolution.
 */
import java.util.Arrays;

public final class AtlasEMStatisticsReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int NX=3, NY=2, NZ=2, NC=2, K=3;
    static final float[][][][] image=new float[NX][NY][NZ][NC];
    static final float[][][][] mem=new float[NX][NY][NZ][K];
    static final float[][][][] field=new float[NX][NY][NZ][NC];
    static final boolean[][][] mask=new boolean[NX][NY][NZ];

    static void init() {
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            mask[x][y][z]=((x+2*y+z)%3)!=0;
            for(int c=0;c<NC;c++) {
                image[x][y][z][c]=(0.25f + (x*17+y*7+z*3+c*11)%29)/7.0f;
                field[x][y][z][c]=0.8f+(x+y+z+c)/10.0f;
            }
            for(int k=0;k<K;k++) mem[x][y][z][k]=(1+x*3+y*2+z+k*5)/17.0f;
        }
        // An empty class checks the source's keep-initial-centroid behavior.
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)mem[x][y][z][2]=0.0f;
    }
    static float[] priors(boolean equal) {
        float[] out=new float[K];
        if(equal) { Arrays.fill(out,1.0f/K); return out; }
        float den=0.0f;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]) {
            den+=1.0f;
            for(int k=0;k<K;k++)out[k]+=mem[x][y][z][k];
        }
        if(den!=0.0f)for(int k=0;k<K;k++)out[k]/=den;
        return out;
    }
    static float[][] centroids() {
        float[][] out={{0.3f,0.4f,0.7f}};
        // actual source shape is channels x classes
        out=new float[NC][K]; out[0][2]=0.7f; out[1][2]=0.8f;
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++) {
            float num=0.0f,den=0.0f;
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]) {
                num+=mem[x][y][z][k]*image[x][y][z][c];
                den+=mem[x][y][z][k];
            }
            if(den!=0.0f)out[c][k]=num/den;
        }
        return out;
    }
    static float[][][] covariance(float[][] means) {
        float[][][] out=new float[K][NC][NC];
        for(int k=0;k<K;k++) {
            float den=0.0f;
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z])
                for(int c1=0;c1<NC;c1++)for(int c2=0;c2<NC;c2++) {
                    out[k][c1][c2]+=(image[x][y][z][c1]-means[c1][k])
                        *(image[x][y][z][c2]-means[c2][k])*mem[x][y][z][k];
                    den+=mem[x][y][z][k];
                }
            for(int c1=0;c1<NC;c1++)for(int c2=0;c2<NC;c2++)if(den!=0.0f)out[k][c1][c2]/=den;
        }
        return out;
    }
    static int bits(float x){return Float.floatToRawIntBits(x);}
    static String vector(float[] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=bits(a[i]);}return s+"]";}
    static String matrix(float[][] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=vector(a[i]);}return s+"]";}
    static String tensor(float[][][] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=matrix(a[i]);}return s+"]";}
    public static void main(String[] args) {
        init(); float[][] mean=centroids();
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"prior\":"+vector(priors(false))+",\"equal_prior\":"+vector(priors(true))
            +",\"centroid\":"+matrix(mean)+",\"covariance\":"+tensor(covariance(mean))+"}");
    }
}
