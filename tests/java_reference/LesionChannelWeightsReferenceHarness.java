/** Selected-method excerpt of LesionToadSegmentation.computeChannelWeights. */
public final class LesionChannelWeightsReferenceHarness {
    private static final int NX=5,NY=4,NZ=3,NC=2,C=6;
    private static final float OUTLIER=2.75f;
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC],field=new float[NX][NY][NZ][NC];
    private static final float[][] centroid={{0.05f,0.10f,0.15f,0.20f,0.0f,0.0f},{0.10f,0.15f,0.20f,0.25f,0.0f,0.0f}};
    private static final float[] lesionCentroid={0.03f,0.04f},blackCentroid={0.02f,0.02f};
    private static final float[][][] lesion=new float[NX][NY][NZ],black=new float[NX][NY][NZ],outlier=new float[NX][NY][NZ];
    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[] compute(String[] types,boolean useField){
        float[] channelWeight=new float[NC];float sum=0.0f;
        for(int c=0;c<NC;c++){
            float num=0.0f,den=0.0f;
            for(int k=0;k<types.length;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
                float m=mem[x][y][z][k],sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                float eps=1e-30f,diff,diffL,diffB,diffO,omem,lmem,bmem;
                if(types[k].equals("object")){
                    diff=sample-centroid[c][k];num+=(m*m+eps)*diff*diff;den+=(m*m+eps);
                }else if(types[k].equals("outlier")){
                    diff=sample-centroid[c][k];diffO=OUTLIER;omem=outlier[x][y][z];
                    num+=(m*m*(1.0f-omem)*(1.0f-omem)+eps)*diff*diff+(m*m*omem*omem)*diffO*diffO;den+=(m*m+eps);
                }else if(types[k].equals("lesion")||types[k].equals("optimized_lesion")||types[k].equals("blackhole")){
                    diff=sample-centroid[c][k];diffL=sample-lesionCentroid[c];diffB=sample-blackCentroid[c];lmem=lesion[x][y][z];bmem=black[x][y][z];
                    num+=(m*m+eps)*diff+(lmem*lmem+eps)*diffL+(bmem*bmem+eps)*diffB;
                    den+=(m*m+eps)+(lmem*lmem+eps)+(bmem*bmem+eps);
                }
            }
            channelWeight[c]=num>0.0f?den/num:1.0e10f;sum+=channelWeight[c];
        }
        for(int c=0;c<NC;c++)channelWeight[c]/=sum;
        return channelWeight;
    }
    private static String vector(float[] a){return "["+bits(a[0])+","+bits(a[1])+"]";}
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
            lesion[x][y][z]=0.1f+((x+2*y+z)%7)*0.05f;
            black[x][y][z]=0.08f+((2*x+y+3*z)%8)*0.04f;
            outlier[x][y][z]=0.05f+((x+3*y+2*z)%7)*0.06f;
        }
        String[] types={"object","outlier","lesion","blackhole","mask","outmask"};
        String[] emptyTypes={"mask","outmask"};
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"diagonal_no_field\":"+vector(compute(types,false))+
            ",\"diagonal_with_field\":"+vector(compute(types,true))+
            ",\"uniform\":["+bits(0.5f)+","+bits(0.5f)+"],\"zero_energy\":"+vector(compute(emptyTypes,false))+"}");
    }
}
