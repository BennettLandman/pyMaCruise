import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected LESION/BLACKHOLE branches of LesionToadSegmentation.computeCentroids(). */
public final class LesionAuxiliaryCentroidsReferenceHarness {
    private static final int NX=5, NY=4, NZ=3, NC=2, C=5;
    private static final float ZERO=1e-10f, SMOOTH=0.35f;
    private static final String SOURCE_SHA="77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0";
    private static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final float[][] INITIAL={{0.15f,0.55f,1.05f,1.65f,2.0f},{0.30f,0.70f,1.15f,1.80f,2.1f}};
    private static final float[][] PRIOR={{0.20f,0.60f,0.95f,1.40f,2.2f},{0.35f,0.75f,1.10f,1.55f,2.3f}};
    private static final float[] LESION_CENTROID={0.85f,0.95f};
    private static final float[] BLACK_CENTROID={0.45f,0.52f};
    private static final float[] LESION_EXP={1.2f,0.8f};
    private static final int[][] LINK={{0,1,1,2,3},{0,0,1,1,2}};
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC];
    private static final float[][][][] field=new float[NX][NY][NZ][NC];
    private static final float[][][] lesion=new float[NX][NY][NZ], black=new float[NX][NY][NZ];
    private static final int[] routeCounts=new int[3], tieCounts=new int[3];

    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[][][] update(String mode, boolean useField){
        float[][] classCentroid=new float[NC][C];
        float[] lesionCentroid=LESION_CENTROID.clone(),blackCentroid=BLACK_CENTROID.clone();
        for(int c=0;c<NC;c++)classCentroid[c]=INITIAL[c].clone();
        PowerTable[] tables=new PowerTable[NC];
        for(int c=0;c<NC;c++)tables[c]=new PowerTable(0.0f,4.0f,0.000001f,LESION_EXP[c]);
        for(int c=0;c<NC;c++){
            float[] num=new float[C],den=new float[C];
            float numL=0.0f,denL=0.0f,numB=0.0f,denB=0.0f;
            for(int k=0;k<C;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                float m=mem[x][y][z][k];
                float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                if(k==1||k==2){
                    float lm=lesion[x][y][z],bm=black[x][y][z],winner=Math.max(m,Math.max(lm,bm));
                    if(winner==m){num[k]+=m*m*sample;den[k]+=m*m;}
                    else if(winner==lm){
                        float residual=sample-lesionCentroid[c],dist=residual*residual,factor=1.0f;
                        if(dist>ZERO)factor=(float)(tables[c].lookup((double)dist,(double)LESION_EXP[c])/dist);
                        numL+=factor*lm*lm*sample;denL+=factor*lm*lm;
                    }else{numB+=bm*bm*sample;denB+=bm*bm;}
                }else{num[k]+=m*m*sample;den[k]+=m*m;}
            }
            for(int k=0;k<C;k++){
                if(mode.equals("none")){
                    if(den[k]>0.0f)classCentroid[c][k]=num[k]/den[k];else classCentroid[c][k]=0.0f;
                }else if(mode.equals("linked")){
                    float ns=num[k],ds=den[k];
                    for(int l=0;l<C;l++)if(l!=k&&LINK[c][k]==LINK[c][l]){ns+=SMOOTH*num[l];ds+=SMOOTH*den[l];}
                    if(ds>0.0f)classCentroid[c][k]=ns/ds;else classCentroid[c][k]=0.0f;
                }else if(mode.equals("damped")){
                    if(den[k]>0.0f)classCentroid[c][k]=(1.0f-SMOOTH)*num[k]/den[k]+SMOOTH*classCentroid[c][k];
                }else if(mode.equals("prior")){
                    float ns=num[k],ds=den[k];ns+=SMOOTH*1e4f*PRIOR[c][k];ds+=SMOOTH*1e4f;
                    if(ds>0.0f)classCentroid[c][k]=ns/ds;else classCentroid[c][k]=0.0f;
                }else{
                    if(den[k]>0.0f)classCentroid[c][k]=num[k]/den[k];else classCentroid[c][k]=0.0f;
                }
            }
            if(denL>0.0f)lesionCentroid[c]=0.2f*(numL/denL)+0.8f*lesionCentroid[c];else lesionCentroid[c]=0.0f;
            if(denB>0.0f)blackCentroid[c]=numB/denB;else blackCentroid[c]=0.0f;
        }
        return new float[][][]{classCentroid,{lesionCentroid},{blackCentroid}};
    }
    private static String matrix(float[][] values){
        StringBuilder b=new StringBuilder("[");
        for(int c=0;c<values.length;c++){if(c>0)b.append(',');b.append('[');for(int k=0;k<values[c].length;k++){if(k>0)b.append(',');b.append(bits(values[c][k]));}b.append(']');}
        return b.append(']').toString();
    }
    private static String vector(float[] values){StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(bits(values[i]));}return b.append(']').toString();}
    private static void updateRoutingEvidence(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int k=seg[x][y][z];if(k!=1&&k!=2)continue;
            float m=mem[x][y][z][k],lm=lesion[x][y][z],bm=black[x][y][z],winner=Math.max(m,Math.max(lm,bm));
            if(m==lm&&m==bm)tieCounts[2]++;
            else if(lm==bm&&lm>m)tieCounts[1]++;
            else if(m==lm&&m>bm)tieCounts[0]++;
            if(winner==m)routeCounts[0]++;else if(winner==lm)routeCounts[1]++;else routeCounts[2]++;
        }
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
            int ls=(x+2*y+z)%9,bs=(2*x+y+3*z)%9;
            lesion[x][y][z]=0.05f+ls*0.06f;
            black[x][y][z]=0.04f+bs*0.07f;
            if((x+y+z)%7==0)black[x][y][z]=lesion[x][y][z];
            if((x+y+z)%13==0)lesion[x][y][z]=black[x][y][z]=0.0f;
        }
        mem[0][0][1][1]=0.8f;lesion[0][0][1]=0.8f;black[0][0][1]=0.2f;
        mem[0][0][2][2]=0.3f;lesion[0][0][2]=0.8f;black[0][0][2]=0.8f;
        mem[1][0][1][2]=0.6f;lesion[1][0][1]=0.6f;black[1][0][1]=0.6f;
        mem[0][1][0][2]=0.3f;lesion[0][1][0]=0.2f;black[0][1][0]=0.9f;
        mem[1][0][0][1]=0.2f;lesion[1][0][0]=0.9f;black[1][0][0]=0.1f;image[1][0][0][0]=0.85f;
        updateRoutingEvidence();
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"powertable_sha256\":\""+POWER_SHA+"\",");
        System.out.print("\"route_counts\":["+routeCounts[0]+","+routeCounts[1]+","+routeCounts[2]+"],\"tie_counts\":["+tieCounts[0]+","+tieCounts[1]+","+tieCounts[2]+"]");
        for(String mode:new String[]{"none","linked","damped","prior","fallback"})for(int fieldMode=0;fieldMode<2;fieldMode++){
            float[][][] o=update(mode,fieldMode==1);String suffix=fieldMode==1?"with_field":"no_field";
            System.out.print(",\""+mode+"_"+suffix+"\":{\"class\":"+matrix(o[0])+",\"lesion\":"+vector(o[1][0])+",\"blackhole\":"+vector(o[2][0])+"}");
        }
        System.out.print("}");
    }
}
