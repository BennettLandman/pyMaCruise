import edu.jhmi.rad.medic.utilities.PowerTable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;

/** Independent selected-method composition for one default exact TOADS slice. */
public final class ToadsDefaultExactIterationSliceReferenceHarness {
    private static final int N=36,K=11,NC=2,BORDER=10;
    private static final float SMOOTH=.043f,FACTOR=.64f,COEFFICIENT=.16f,PRIOR_SCALE=.05f,ZERO=1e-10f,INF=1e10f;
    private static final float[] VARIMG={.61f,.34f},MODALITY={.73f,.46f};
    private static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] VAR_WEIGHT={
        {.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},
        {1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}
    };
    private static final float[][] EXP={{1,.8f,.8f,1,1,1,1,1,1,1,1},{1,.8f,.8f,1,1,1,1,1,1,1,1}};
    private static final int[] CATEGORY={-1,0,0,1,1,1,1,1,2,2,2};
    private static final float[][][][] image=new float[N][N][N][NC],field=new float[N][N][N][NC];
    private static final int[][][] segmentation=new int[N][N][N];
    private static final float[][][][] mem=new float[N][N][N][K],shapes=new float[N][N][N][K];
    private static final int[][][][] closest=new int[N][N][N][4];
    private static final float[][] center=new float[NC][K],std=new float[NC][K],similarity=new float[K][K];
    private static final float[][][][] energy=new float[N][N][N][K];
    private static final float[][][] denominator=new float[N][N][N];
    private static final PowerTable[][] powers=new PowerTable[3][NC];

    private static int label(int x,int y,int z){return Math.floorMod(x+2*y+3*z+(y*z)%3+(x*y)%2,K);}
    private static float imageValue(int channel,int x,int y,int z,int k){
        int delta=channel==0?(x*13+y*7+z*3)%17-8:(x*5+y*11+z*7)%19-9;
        float base=channel==0?(k+1)/12.0f:(11-k)/12.0f;
        return base+(channel==0?delta/1700.0f:delta/1900.0f);
    }
    private static float fieldValue(int x,int y,int z){return 1.0f+((x+3*y+5*z)%9)/100.0f;}
    private static float shapeValue(int x,int y,int z,int k){return .20f+((7*x+5*y+3*z+2*k)%70)/100.0f;}
    private static int[] closestValue(int x,int y,int z){return new int[]{(x+y+z)%K,(x+2*y+3*z+1)%K,(x+3*y+2*z+4)%K,(2*x+y+z+7)%K};}
    private static void initialize(){
        float[][] old=new float[NC][K];
        for(int k=0;k<K;k++){old[0][k]=(k+1)/12.0f;old[1][k]=(11-k)/12.0f;}
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int current=label(x,y,z);segmentation[x][y][z]=current;
            for(int c=0;c<NC;c++){image[x][y][z][c]=imageValue(c,x,y,z,current);field[x][y][z][c]=fieldValue(x,y,z);}
            for(int k=0;k<K;k++)shapes[x][y][z][k]=shapeValue(x,y,z,k);
            closest[x][y][z]=closestValue(x,y,z);
        }
        for(int x=BORDER;x<N-BORDER-1;x++)for(int y=BORDER;y<N-BORDER-1;y++)for(int z=BORDER;z<N-BORDER-1;z++){
            float den=0.0f;
            for(int k=0;k<K;k++){
                float d=0.0f;
                for(int c=0;c<NC;c++){float difference=image[x][y][z][c]-old[c][k];d+=difference*difference;}
                float v=d>1e-20f?1.0f/d:1e20f;mem[x][y][z][k]=v;den+=v;
            }
            if(den>0.0f)for(int k=0;k<K;k++)mem[x][y][z][k]/=den;
        }
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++){
            float numerator=0.0f,den=0.0f;
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
                float m=mem[x][y][z][k],weight=m*m,corrected=field[x][y][z][c]*image[x][y][z][c];
                float w=weight;
                if(k>0){float residual=corrected-old[c][k],distance=residual*residual;float exponent=EXP[c][k];
                    float factor=1.0f;if(distance>1e-10f)factor=(float)(powersFor(c,k).lookup((double)distance,(double)exponent)/distance);
                    w=factor*m;w*=m;
                }
                numerator+=w*corrected;den+=w;
            }
            center[c][k]=den>0.0f?(1.0f-.1f)*numerator/den+.1f*old[c][k]:old[c][k];
        }
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++){
            float num=0.0f,den=0.0f;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)if(segmentation[x][y][z]==k){
                float m=mem[x][y][z][k],weight=m*m+1e-30f;
                float corrected=field[x][y][z][c]*image[x][y][z][c],diff=corrected-center[c][k];
                num+=weight*diff*diff;den+=weight;
            }
            std[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
        for(int t=0;t<3;t++)for(int c=0;c<NC;c++){
            float exponent=t==0?.8f:1.0f;powers[t][c]=new PowerTable(0,4,.000001f,exponent);
        }
    }
    private static PowerTable powersFor(int c,int k){return powers[CATEGORY[k]][c];}
    private static float shapePriorSimilarity(int k,int m){
        if(k==m)return 0.0f;float distance=0.0f;
        for(int c=0;c<NC;c++){
            float delta=center[c][k]-center[c][m],sum=std[c][k]+std[c][m];
            float term=VARIMG[c]*delta;term*=delta;term/=sum;term/=sum;distance+=term;
        }
        float scale=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;
        return 1.0f/(1.0f+scale*distance);
    }
    private static boolean direct(int x,int y,int z,int k){return closest[x][y][z][0]==k||closest[x][y][z][1]==k;}
    private static boolean neighbor(int x,int y,int z,int k){return closest[x][y][z][2]==k||closest[x][y][z][3]==k;}
    private static float relation(int x,int y,int z,int k){float s=shapes[x][y][z][k];if(direct(x,y,z,k))return 1.0f;if(neighbor(x,y,z,k))return Math.max(.25f,Math.min(2.0f*FACTOR*s,1.0f));return Math.min(1.0f,2.0f*FACTOR*s);}
    private static void update(){
        for(int k=0;k<K;k++)for(int m=0;m<K;m++)similarity[k][m]=shapePriorSimilarity(k,m);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<K;k++){
            if(CATEGORY[k]<0){for(int c=0;c<NC;c++){float w=MODALITY[c]*VAR_WEIGHT[c][k];w*=VARIMG[c];float d=image[x][y][z][c];energy[x][y][z][k]+=w*(d==0.0f?0.0f:1.0f);}}
            else for(int c=0;c<NC;c++){
                float corrected=field[x][y][z][c]*image[x][y][z][c],residual=corrected-center[c][k],d=residual*residual;
                double term=powersFor(c,k).lookup((double)d,(double)EXP[c][k]);float w=MODALITY[c]*VAR_WEIGHT[c][k];w*=VARIMG[c];
                energy[x][y][z][k]=(float)(energy[x][y][z][k]+w*term);
            }
        }
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            float den=0.0f;float[] objectiveAt=new float[K];
            for(int k=0;k<K;k++){
                boolean d=direct(x,y,z,k),n=neighbor(x,y,z,k);if(!d&&!n&&shapes[x][y][z][k]==0.0f){mem[x][y][z][k]=0.0f;continue;}
                float num=energy[x][y][z][k],obj=num;
                float neighborSum=0.0f;
                for(int[] offset:offsets)for(int m=0;m<K;m++)if(m!=k){float value=mem[x+offset[0]][y+offset[1]][z+offset[2]][m];neighborSum+=value*value;}
                float increment=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;num+=increment;obj+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                float inversePrior=0.0f;
                for(int m=0;m<K;m++)if(m!=k){float value=ATLAS_PRIOR[m]*similarity[k][m];value*=shapes[x][y][z][m];value*=shapes[x][y][z][m];inversePrior+=value;}
                num+=COEFFICIENT*inversePrior;obj+=.5f*COEFFICIENT*inversePrior;
                float reciprocal=num>ZERO?1.0f/num:INF;mem[x][y][z][k]=relation(x,y,z,k)*reciprocal;den+=mem[x][y][z][k];objectiveAt[k]=obj;
            }
            denominator[x][y][z]=den;
            if(den>0.0f)for(int k=0;k<K;k++){mem[x][y][z][k]/=den;float contribution=objectiveAt[k]*mem[x][y][z][k];objective[0]+=contribution*mem[x][y][z][k];}
        }
    }
    private static final float[] objective={0.0f};
    private static String hash(float[][][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hash(float[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][] x:a)for(float[] y:x)for(float v:y){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hash(float[][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[] x:a)for(float v:x){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hex(byte[] a){StringBuilder s=new StringBuilder();for(byte b:a)s.append(String.format("%02x",b));return s.toString();}
    private static String bits(float[][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(float[] row:a)for(float v:row){if(!first)s.append(',');first=false;s.append(Float.floatToIntBits(v));}return s.append(']').toString();}
    public static void main(String[] args)throws Exception{
        for(int c=0;c<NC;c++)for(int t=0;t<3;t++)powers[t][c]=new PowerTable(0,4,.000001f,t==0?.8f:1.0f);
        initialize();update();
        System.out.println("{\"source\":\"VandyCRUISE MedicAlgorithmToads08 exact default startup -> computeExactMemberships selected-method composition\",\"segmentation_source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"preprocess_source_sha256\":\"31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7\",\"power_table_source_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"shape_xyz_channels\":[36,36,36,2],\"classes\":11,\"explicit_image_variance_weights\":[0.61,0.34],\"smoothing\":0.043,\"prior_factor\":0.64,\"prior_coefficient\":0.16,\"prior_scale\":0.05,\"centroids_bits\":"+bits(center)+",\"standard_deviation_bits\":"+bits(std)+",\"similarity_sha256\":\""+hash(similarity)+"\",\"class_energy_sha256\":\""+hash(energy)+"\",\"denominator_sha256\":\""+hash(denominator)+"\",\"membership_bits_sha256\":\""+hash(mem)+"\",\"objective_bits\":"+Float.floatToIntBits(objective[0])+"}");
    }
}
