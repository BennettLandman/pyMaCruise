/** Selected-method excerpt of LesionToadSegmentation.computeVariances. */
public final class LesionVariancesReferenceHarness {
    private static final int NX=5,NY=4,NZ=3,NC=2,C=4;
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC];
    private static final float[][][][] field=new float[NX][NY][NZ][NC];
    private static final float[][] centroid={{0.15f,0.55f,1.05f,1.65f},{0.30f,0.70f,1.15f,1.80f}};
    private static final float[][] lesionCentroid={{0.85f},{0.95f}};
    private static final float[][] blackCentroid={{0.45f},{0.52f}};
    private static final float[][][] lesion=new float[NX][NY][NZ],black=new float[NX][NY][NZ];
    private static final String[] type={"object","lesion","blackhole","object"};

    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static final class Result {float[][] cls=new float[NC][C];float[] les=new float[NC],blk=new float[NC];}
    private static Result update(boolean useField){
        Result out=new Result();
        for(int c=0;c<NC;c++){
            float lnum=0.0f,lden=0.0f,bnum=0.0f,bden=0.0f;
            for(int k=0;k<C;k++){
                float num=0.0f,den=0.0f;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                    float m=mem[x][y][z][k],lm=lesion[x][y][z],bm=black[x][y][z];
                    boolean special=((type[k].equals("lesion")||type[k].equals("optimized_lesion"))&&m<lm)||(type[k].equals("blackhole")&&m<bm);
                    float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                    if(special){
                        float selected=lm;
                        if(lm>=bm){float diff=sample-lesionCentroid[c][0];lnum+=(selected*selected+1e-30f)*diff*diff;lden+=selected*selected+1e-30f;}
                        else {float diff=sample-blackCentroid[c][0];bnum+=(selected*selected+1e-30f)*diff*diff;bden+=selected*selected+1e-30f;}
                    }else{
                        float diff=sample-centroid[c][k];num+=(m*m+1e-30f)*diff*diff;den+=m*m+1e-30f;
                    }
                }
                out.cls[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
            }
            out.les[c]=lden>0.0f?(float)Math.sqrt(lnum/lden):0.0f;
            out.blk[c]=bden>0.0f?(float)Math.sqrt(bnum/bden):0.0f;
        }
        return out;
    }
    private static String matrix(float[][] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('[');for(int j=0;j<a[i].length;j++){if(j>0)b.append(',');b.append(bits(a[i][j]));}b.append(']');}return b.append(']').toString();}
    private static String vector(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();}
    private static String json(Result r){return "{\"class_bits\":"+matrix(r.cls)+",\"lesion_bits\":"+vector(r.les)+",\"black_hole_bits\":"+vector(r.blk)+"}";}
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            seg[x][y][z]=(short)((x+2*y+z)%3);
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
            float lm=0.05f+((x+2*y+z)%9)*0.06f;
            float bm=0.04f+((2*x+y+3*z)%9)*0.07f;
            if((x+y+z)%7==0)bm=lm;
            if((x+y+z)%13==0){lm=0.0f;bm=0.0f;}
            lesion[x][y][z]=lm;black[x][y][z]=bm;
        }
        Result raw=update(false),corrected=update(true);
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"no_field\":"+json(raw)+",\"with_field\":"+json(corrected)+"}");
    }
}
