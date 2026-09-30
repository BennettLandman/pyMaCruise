/** Selected generic branch of LesionToadSegmentation.computeCentroids(). */
public final class LesionBasicCentroidsReferenceHarness {
    private static final int NX=5, NY=4, NZ=3, NC=2, C=4;
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC];
    private static final float[][][][] field=new float[NX][NY][NZ][NC];
    private static final float[][] INITIAL={{0.15f,0.55f,1.05f,1.65f},{0.30f,0.70f,1.15f,1.80f}};
    private static final float[][] PRIOR={{0.20f,0.60f,0.95f,1.40f},{0.35f,0.75f,1.10f,1.55f}};
    private static final int[][] LINK={{0,1,1,2},{0,0,1,1}};
    private static final float SMOOTH=0.35f;

    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[][] update(String mode, boolean useField){
        float[][] centroid=new float[NC][C];
        for(int c=0;c<NC;c++)centroid[c]=INITIAL[c].clone();
        for(int c=0;c<NC;c++){
            float[] num=new float[C], den=new float[C];
            for(int k=0;k<C;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                float m=mem[x][y][z][k];
                if(k==1){
                    if(useField)num[k]+=m*m*(field[x][y][z][c]*image[x][y][z][c])*(field[x][y][z][c]*image[x][y][z][c]);
                    else num[k]+=m*m*(image[x][y][z][c])*(image[x][y][z][c]);
                }else{
                    if(useField)num[k]+=m*m*(field[x][y][z][c]*image[x][y][z][c]);
                    else num[k]+=m*m*(image[x][y][z][c]);
                }
                den[k]+=m*m;
            }
            for(int k=0;k<C;k++){
                boolean rayleigh=k==1;
                if(mode.equals("none")){
                    if(den[k]>0.0f){
                        if(rayleigh)centroid[c][k]=(float)Math.sqrt(num[k]/(2.0f*den[k]));
                        else centroid[c][k]=num[k]/den[k];
                    }else centroid[c][k]=0.0f;
                }else if(mode.equals("linked")){
                    float nsum=num[k],dsum=den[k];
                    for(int l=0;l<C;l++)if(l!=k&&LINK[c][k]==LINK[c][l]){
                        nsum+=SMOOTH*num[l];dsum+=SMOOTH*den[l];
                    }
                    if(dsum>0.0f){
                        if(rayleigh)centroid[c][k]=(float)Math.sqrt(nsum/(2.0f*dsum));
                        else centroid[c][k]=nsum/dsum;
                    }else centroid[c][k]=0.0f;
                }else if(mode.equals("damped")){
                    if(den[k]>0.0f){
                        if(rayleigh)centroid[c][k]=(1.0f-SMOOTH)*(float)Math.sqrt(num[k]/(2.0f*den[k]))+SMOOTH*centroid[c][k];
                        else centroid[c][k]=(1.0f-SMOOTH)*num[k]/den[k]+SMOOTH*centroid[c][k];
                    }
                }else if(mode.equals("prior")){
                    float nsum=num[k],dsum=den[k];
                    nsum+=SMOOTH*1e4f*PRIOR[c][k];dsum+=SMOOTH*1e4f;
                    if(dsum>0.0f)centroid[c][k]=nsum/dsum;
                    else centroid[c][k]=0.0f;
                }else{
                    if(den[k]>0.0f){
                        if(rayleigh)centroid[c][k]=(float)Math.sqrt(num[k]/(2.0f*den[k]));
                        else centroid[c][k]=num[k]/den[k];
                    }else centroid[c][k]=0.0f;
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
            seg[x][y][z]=(short)((x+2*y+z)%3);
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
        }
        System.out.print("{\"source_sha256\":\"77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0\"");
        for(String mode:new String[]{"none","linked","damped","prior","fallback"}){
            System.out.print(",\""+mode+"_no_field\":"+bitsJson(update(mode,false)));
            System.out.print(",\""+mode+"_with_field\":"+bitsJson(update(mode,true)));
        }
        System.out.print("}");
    }
}
