import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** P18czo: run the selected synthetic two-channel continuation with supplied fields. */
public final class ToadsDefault11DynamicMultichannelMembershipReferenceHarness extends ToadsDefault11ContinueThinningReferenceHarness {
    private static final int NC=2;
    private static final float SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f, PRIOR_SCALE=.05f;
    private static final float[] MODALITY={.73f,.46f}, IMAGE_VARIANCE={.61f,.34f};
    private static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] VARIANCE_WEIGHT={
        {.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},
        {1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}
    };
    private static final float[][] EXPONENT=new float[NC][11];
    private static final int[] CATEGORY={-1,0,0,1,1,1,1,1,2,2,2};
    private final float[][][][] images2=new float[N][N][N][NC], fields2=new float[N][N][N][NC];
    private final float[][] centers=new float[NC][C], deviations=new float[NC][C];
    private final float[][] similarity=new float[C][C];
    private final float[][][][] energy=new float[N][N][N][C];
    private final float[][][] denominator=new float[N][N][N];
    private final PowerTable[][] powers=new PowerTable[3][NC];
    private float objective=0.0f;

    private static float secondImage(int x,int y,int z) {
        int label=(x*3+y*5+z*7+(x*y)%C)%C;
        for(int k=0;k<C;k++)if(x==1+k&&y==1+(k*3)%12&&z==1+(k*5)%12)label=k;
        int delta=(x*5+y*11+z*7)%19-9;
        return (11-label)/12.0f+delta/1900.0f;
    }
    private static int[] closest(int x,int y,int z) {
        return new int[]{(x+y+z)%C,(x+2*y+3*z+1)%C,(x+3*y+2*z+4)%C,(2*x+y+z+7)%C};
    }
    private void loadTwoChannels(Path fieldPath) throws IOException {
        byte[] bytes=Files.readAllBytes(fieldPath);
        int count=N*N*N*NC;
        if(bytes.length!=count*Float.BYTES)throw new IllegalArgumentException("expected 14x14x14x2 little-endian float32 fields");
        ByteBuffer input=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            float f0=input.getFloat(),f1=input.getFloat();
            field[x][y][z][0]=f0;
            fields2[x][y][z][0]=f0;
            fields2[x][y][z][1]=f1;
            images2[x][y][z][0]=images[x][y][z][0];
            images2[x][y][z][1]=secondImage(x,y,z);
        }
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)EXPONENT[c][k]=CATEGORY[k]<0?1.0f:(CATEGORY[k]==0?.8f:1.0f);
    }
    private void refreshMultichannelStatistics() {
        PowerTable[] tables={new PowerTable(0,4,.000001f,.8f),new PowerTable(0,4,.000001f,1.0f),new PowerTable(0,4,.000001f,1.0f)};
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++) {
            float old=c==0?oldCentroids[k]:(11-k)/12.0f, num=0.0f, den=0.0f;
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c], mem=mems[x][y][z][k], weight=mem*mem;
                if(k>0) {
                    float residual=corrected-old,distance=residual*residual,factor=1.0f;
                    if(distance>1e-10f)factor=(float)(tables[CATEGORY[k]].lookup((double)distance,(double)EXPONENT[c][k])/distance);
                    weight=factor*mem;weight*=mem;
                }
                num+=weight*corrected;den+=weight;
            }
            centers[c][k]=den>0.0f?.9f*num/den+.1f*old:old;
        }
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++) {
            float num=0.0f,den=0.0f;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)if(seg[x][y][z]==k) {
                float mem=mems[x][y][z][k],weight=mem*mem+1e-30f;
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c],diff=corrected-centers[c][k];
                num+=weight*diff*diff;den+=weight;
            }
            deviations[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
    }
    private void updateMemberships() {
        for(int type=0;type<3;type++)for(int c=0;c<NC;c++)powers[type][c]=new PowerTable(0,4,.000001f,type==0?.8f:1.0f);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<C;k++) {
            if(CATEGORY[k]<0) {
                for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=IMAGE_VARIANCE[c];float v=images2[x][y][z][c];if(v!=0.0f)energy[x][y][z][k]+=w;}
            } else for(int c=0;c<NC;c++) {
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c],residual=corrected-centers[c][k],d=residual*residual;
                double term=powers[CATEGORY[k]][c].lookup((double)d,(double)EXPONENT[c][k]);
                float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=IMAGE_VARIANCE[c];
                energy[x][y][z][k]=(float)(energy[x][y][z][k]+w*term);
            }
        }
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int[] labels=closest(x,y,z);float den=0.0f;float[] objectiveAt=new float[C];
            for(int k=0;k<C;k++) {
                boolean direct=labels[0]==k||labels[1]==k,neighbor=labels[2]==k||labels[3]==k;
                if(!direct&&!neighbor&&shape[x][y][z][k]==0.0f){energy[x][y][z][k]=0.0f;mems[x][y][z][k]=0.0f;continue;}
                float num=energy[x][y][z][k],obj=num,neighborSum=0.0f;
                for(int[] d:offsets)for(int m=0;m<C;m++)if(m!=k){float v=mems[x+d[0]][y+d[1]][z+d[2]][m];neighborSum+=v*v;}
                float increment=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;num+=increment;obj+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                float inverse=0.0f;
                for(int m=0;m<C;m++)if(m!=k){float v=ATLAS_PRIOR[m]*similarity[k][m];v*=shape[x][y][z][m];v*=shape[x][y][z][m];inverse+=v;}
                num+=COEFFICIENT*inverse;obj+=.5f*COEFFICIENT*inverse;
                float factor=shape[x][y][z][k];
                if(direct)factor=1.0f;else if(neighbor)factor=Math.max(.25f,Math.min(2.0f*FACTOR*factor,1.0f));else factor=Math.min(1.0f,2.0f*FACTOR*factor);
                float value=factor*(num>1e-10f?1.0f/num:1e10f);mems[x][y][z][k]=value;den+=value;objectiveAt[k]=obj;
            }
            denominator[x][y][z]=den;
            if(den>0.0f)for(int k=0;k<C;k++){mems[x][y][z][k]/=den;float contribution=objectiveAt[k]*mems[x][y][z][k];objective+=contribution*mems[x][y][z][k];}
        }
    }
    private static String hash(float[][][][] a)throws Exception {MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hash(float[][][] a)throws Exception {MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][] x:a)for(float[] y:x)for(float v:y){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hex(byte[] a){StringBuilder s=new StringBuilder();for(byte b:a)s.append(String.format("%02x",b));return s.toString();}
    private static String labels(byte[][][] a){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(v<10?(char)('0'+v):(char)('a'+v-10));return b.toString();}
    private static String bits(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();}
    public static void main(String[] args) {
        if(args.length!=1)throw new IllegalArgumentException("usage: harness two-channel-fields.bin");
        try {
            ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
            ToadsDefault11DynamicMultichannelMembershipReferenceHarness h=new ToadsDefault11DynamicMultichannelMembershipReferenceHarness();
            h.initialize();h.grow();h.refreshRelationAndInitializeThinning();h.thin();h.loadTwoChannels(Path.of(args[0]));h.refreshMultichannelStatistics();h.updateMemberships();
            System.out.print("{\"reference_scope\":\"selected Java formula harness; source-connected synthetic grow/thin state, supplied two-channel fields, fixed modality/image/class weights and zero shape maps\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"power_table_source_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"segmentation\":\""+labels(h.seg)+"\",\"centroids_channel0_bits\":"+bits(h.centers[0])+",\"centroids_channel1_bits\":"+bits(h.centers[1])+",\"deviations_channel0_bits\":"+bits(h.deviations[0])+",\"deviations_channel1_bits\":"+bits(h.deviations[1])+",\"membership_sha256\":\""+hash(h.mems)+"\",\"denominator_sha256\":\""+hash(h.denominator)+"\",\"class_energy_sha256\":\""+hash(h.energy)+"\",\"objective_bits\":"+Float.floatToIntBits(h.objective)+"}");
        } catch(Exception e) { throw new RuntimeException(e); }
    }
}
