/** Selected-method excerpt of LesionToadSegmentation.computeNoiseVariance. */
public final class LesionNoiseVarianceReferenceHarness {
    private static final int NX=5,NY=4,NZ=3,NC=2,C=4;
    private static final short[][][] seg=new short[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C];
    private static final float[][][][] image=new float[NX][NY][NZ][NC],field=new float[NX][NY][NZ][NC];
    private static final float[][] centroid={{0.15f,0.55f,1.05f,1.65f},{0.30f,0.70f,1.15f,1.80f}};
    private static int bits(float v){return Float.floatToRawIntBits(v);}
    private static float[][] update(boolean useField,boolean empty){
        float[][] stddev=new float[NC][C];
        for(int c=0;c<NC;c++){
            float num=0.0f,den=0.0f;
            for(int k=1;k<C;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(!empty&&seg[x][y][z]==k){
                float m=mem[x][y][z][k];float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                float diff=sample-centroid[c][k];num+=(m*m+1e-30f)*diff*diff;den+=(m*m+1e-30f);
            }
            if(den>0.0f)for(int k=0;k<C;k++)stddev[c][k]=(float)Math.sqrt(num/den);
            else for(int k=0;k<C;k++)stddev[c][k]=0.0f;
        }
        return stddev;
    }
    private static String matrix(float[][] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('[');for(int j=0;j<a[i].length;j++){if(j>0)b.append(',');b.append(bits(a[i][j]));}b.append(']');}return b.append(']').toString();}
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            seg[x][y][z]=(short)((x+2*y+z)%3);
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.1f+((x*3+y*5+z*7+k*2)%9)*0.05f;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.2f+((x*11+y*7+z*3+c*5)%23)*0.1f;
                field[x][y][z][c]=0.8f+((x*2+y*3+z*4+c)%5)*0.1f;
            }
        }
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"no_field\":"+matrix(update(false,false))+
            ",\"with_field\":"+matrix(update(true,false)) + ",\"empty_support\":"+matrix(update(false,true))+"}");
    }
}
