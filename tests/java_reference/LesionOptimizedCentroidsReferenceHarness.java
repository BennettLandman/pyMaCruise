import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected optimized CSF/GM/WM branches of LesionToadSegmentation.computeCentroids(). */
public final class LesionOptimizedCentroidsReferenceHarness {
    private static final int NX=5, NY=4, NZ=3, NC=2, C=5;
    private static final float ZERO=1e-10f, SMOOTH=0.35f;
    private static final String SOURCE_SHA="77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0";
    private static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final float[][] INITIAL={{0.15f,0.55f,1.05f,1.65f,2.0f},{0.30f,0.70f,1.15f,1.80f,2.1f}};
    private static final float[][] PRIOR={{0.20f,0.60f,0.95f,1.40f,2.2f},{0.35f,0.75f,1.10f,1.55f,2.3f}};
    private static final float[][] EXP={{0.0f,1.35f,0.8f,1.2f,0.0f},{0.0f,1.1f,1.25f,0.75f,0.0f}};
    private static final int[][] LINK={{0,1,1,2,3},{0,0,1,1,2}};
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC];
    private static final float[][][][] field=new float[NX][NY][NZ][NC];

    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[][] update(String mode, boolean useField){
        float[][] centroid=new float[NC][C];
        for(int c=0;c<NC;c++)centroid[c]=INITIAL[c].clone();
        PowerTable[][] tables=new PowerTable[NC][C];
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)
            tables[c][k]=new PowerTable(0.0f,4.0f,0.000001f,EXP[c][k]);
        for(int c=0;c<NC;c++){
            float[] num=new float[C],den=new float[C];
            for(int k=0;k<C;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                boolean optimized=k==2||k==3||(k==1&&c==0);
                float factor=1.0f;
                if(optimized){
                    float residual=sample-centroid[c][k];
                    float dist=residual*residual;
                    if(dist>ZERO)factor=(float)(tables[c][k].lookup((double)dist,(double)EXP[c][k])/dist);
                }
                float m=mem[x][y][z][k];
                num[k]+=factor*m*m*sample;
                den[k]+=factor*m*m;
            }
            for(int k=0;k<C;k++){
                if(mode.equals("none")){
                    if(den[k]>0.0f)centroid[c][k]=num[k]/den[k];else centroid[c][k]=0.0f;
                }else if(mode.equals("linked")){
                    float nsum=num[k],dsum=den[k];
                    for(int l=0;l<C;l++)if(l!=k&&LINK[c][k]==LINK[c][l]){nsum+=SMOOTH*num[l];dsum+=SMOOTH*den[l];}
                    if(dsum>0.0f)centroid[c][k]=nsum/dsum;else centroid[c][k]=0.0f;
                }else if(mode.equals("damped")){
                    if(den[k]>0.0f)centroid[c][k]=(1.0f-SMOOTH)*num[k]/den[k]+SMOOTH*centroid[c][k];
                }else if(mode.equals("prior")){
                    float nsum=num[k],dsum=den[k];
                    nsum+=SMOOTH*1e4f*PRIOR[c][k];dsum+=SMOOTH*1e4f;
                    if(dsum>0.0f)centroid[c][k]=nsum/dsum;else centroid[c][k]=0.0f;
                }else{
                    if(den[k]>0.0f)centroid[c][k]=num[k]/den[k];else centroid[c][k]=0.0f;
                }
            }
        }
        return centroid;
    }
    private static String bitsJson(float[][] values){
        StringBuilder b=new StringBuilder("[");
        for(int c=0;c<NC;c++){if(c>0)b.append(',');b.append('[');for(int k=0;k<C;k++){if(k>0)b.append(',');b.append(bits(values[c][k]));}b.append(']');}
        return b.append(']').toString();
    }
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            seg[x][y][z]=(short)((x+2*y+z)%4);
            for(int k=0;k<4;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            mem[x][y][z][4]=0.1f+((x*5+y*7+z*3+2)%11)*0.04f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
        }
        image[1][0][0][0]=0.55f;
        image[0][0][2][0]=1.050001f;
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"powertable_sha256\":\""+POWER_SHA+"\"");
        for(String mode:new String[]{"none","linked","damped","prior","fallback"}){
            System.out.print(",\""+mode+"_no_field\":"+bitsJson(update(mode,false)));
            System.out.print(",\""+mode+"_with_field\":"+bitsJson(update(mode,true)));
        }
        System.out.print("}");
    }
}
