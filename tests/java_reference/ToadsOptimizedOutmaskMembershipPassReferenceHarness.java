import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships excerpt: optimized tissues plus OUTMASK. */
public final class ToadsOptimizedOutmaskMembershipPassReferenceHarness {
    private static final String TOADS_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWER_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int NX=9,NY=8,NZ=7,NC=2,K=6,OUT=1;
    private static final String[] TYPES={"OPTIMIZED_GM","OUTMASK","MASK","OBJECT","OPTIMIZED_CSF","OPTIMIZED_WM"};
    private static final float[][] CENTROID={{.30f,.34f,.38f,.42f,.46f,.50f},{.26f,.30f,.34f,.38f,.42f,.46f}};
    private static final float[][] VARIANCE={{.82f,1.17f,.69f,.95f,1.23f,.76f},{1.12f,.77f,.91f,1.04f,.68f,1.31f}};
    private static final float[] MODALITY={.73f,.46f},IMAGE_VARIANCE={.61f,.34f};
    private static final float[] EXP_CSF={.78f,.91f},EXP_GM={1.16f,1.29f},EXP_WM={1.55f,1.73f};
    private static final float[] SMOOTH_PRIOR={.65f,.91f,1.13f,.72f,.87f,1.05f};
    private static final float[] ATLAS_PRIOR={.15f,.21f,.12f,.19f,.16f,.17f};
    private static final float[][] SIM={{0,.17f,.24f,.11f,.19f,.14f},{.13f,0,.19f,.20f,.18f,.15f},{.22f,.16f,0,.12f,.20f,.17f},{.16f,.18f,.13f,0,.12f,.16f},{.18f,.14f,.19f,.17f,0,.11f},{.19f,.15f,.17f,.14f,.13f,0}};
    private static final float OUTLIER_VALUE=1.37f,SMOOTH=.043f,FACTOR=.64f,COEFFICIENT=.16f;
    private static final float ZERO=1e-10f,INF=1e10f;
    private static float[][][][] images,field,bias,initial,shapes;
    private static float[][][] initialOutlier;
    private static short[][][][] labels;
    private static PowerTable[][] power;

    private static byte[] unzip(String path)throws Exception{ByteArrayOutputStream o=new ByteArrayOutputStream();try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] a=new byte[8192];int n;while((n=in.read(a))>=0)o.write(a,0,n);}return o.toByteArray();}
    private static float[][][][] readFloat4(String path,int nc)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=nc)throw new IllegalArgumentException("unexpected float NIfTI");int o=(int)b.getFloat(108);float[][][][] a=new float[NX][NY][NZ][nc];for(int c=0;c<nc;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getFloat(o+4*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static float[][][] readFloat3(String path)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=3||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ)throw new IllegalArgumentException("unexpected float3 NIfTI");int o=(int)b.getFloat(108);float[][][] a=new float[NX][NY][NZ];for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z]=b.getFloat(o+4*((z*NY+y)*NX+x));return a;}
    private static short[][][][] readShort4(String path)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=4)throw new IllegalArgumentException("unexpected int16 NIfTI");int o=(int)b.getFloat(108);short[][][][] a=new short[NX][NY][NZ][4];for(int c=0;c<4;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getShort(o+2*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static float exponent(int k,int c){if(TYPES[k].equals("OPTIMIZED_CSF"))return EXP_CSF[c];if(TYPES[k].equals("OPTIMIZED_GM"))return EXP_GM[c];return EXP_WM[c];}
    private static void initTables(){power=new PowerTable[K][NC];for(int k=0;k<K;k++)if(TYPES[k].startsWith("OPTIMIZED_"))for(int c=0;c<NC;c++)power[k][c]=new PowerTable(0.0f,4.0f,0.000001f,exponent(k,c));}
    private static float[][][][] dataEnergy(boolean useField,boolean useBias,float[][][] aux){float[][][][] e=new float[NX][NY][NZ][K];for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<K;k++){
        float n=0.0f,nO=0.0f;
        if(TYPES[k].equals("OBJECT"))for(int c=0;c<NC;c++){float value=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];if(useBias)value-=bias[x][y][z][c];float d=value-CENTROID[c][k],w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];n+=w*d*d;}
        else if(TYPES[k].equals("MASK")||TYPES[k].equals("OUTMASK"))for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];n+=w*(images[x][y][z][c]==0.0f?0.0f:1.0f);if(TYPES[k].equals("OUTMASK")){nO+=w*OUTLIER_VALUE;aux[x][y][z]=nO;}}
        else for(int c=0;c<NC;c++){float value=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];if(useBias)value-=bias[x][y][z][c];float d=value-CENTROID[c][k];float dist=d*d;float w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];n=(float)(n+(double)w*power[k][c].lookup((double)dist,(double)exponent(k,c)));}
        e[x][y][z][k]=n;
    }return e;}
    private static float factor(float shape,int x,int y,int z,int k){boolean direct=labels[x][y][z][0]==k||labels[x][y][z][1]==k,neighbor=labels[x][y][z][2]==k||labels[x][y][z][3]==k;if(direct)return 1.0f;if(neighbor)return Math.max(.25f,Math.min(2.0f*FACTOR*shape,1.0f));return Math.min(1.0f,2.0f*FACTOR*shape);}
    private static float inversePrior(int x,int y,int z,int k){float v=0.0f;for(int m=0;m<K;m++)if(m!=k)v+=ATLAS_PRIOR[m]*SIM[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];return v;}
    private static float[][][] lastOutlier;private static float[][][] lastDen;
    private static float[][][][] run(boolean useField,boolean useBias,boolean useSmooth,boolean usePrior,float[][][][] energy,float[][][] aux,float[] objective){float[][][][] mem=new float[NX][NY][NZ][K];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++)mem[x][y][z][k]=initial[x][y][z][k];float[][][] out=new float[NX][NY][NZ];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)out[x][y][z]=initialOutlier[x][y][z];float[][][] denMap=new float[NX][NY][NZ];float obj=0.0f;int[][] ds={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){float den=0.0f;float[] objAt=new float[K];for(int k=0;k<K;k++){float num=energy[x][y][z][k],numO=aux[x][y][z],objectiveValue=num,oldOut=out[x][y][z];if(k==OUT)objectiveValue=(1.0f-oldOut)*num+oldOut*numO;
            if(useSmooth){float ngb=0.0f;for(int[] d:ds)for(int m=0;m<K;m++)if(m!=k){float v=mem[x+d[0]][y+d[1]][z+d[2]][m];ngb+=v*v;}float inc=SMOOTH_PRIOR[k]*SMOOTH*ngb;num+=inc;objectiveValue+=.5f*SMOOTH_PRIOR[k]*SMOOTH*ngb;/* OUTMASK auxiliary smoothing skipped by source dangling else. */}
            float inv=usePrior?inversePrior(x,y,z,k):0.0f;if(usePrior){num+=COEFFICIENT*inv;if(k==OUT)numO+=COEFFICIENT*inv;objectiveValue+=.5f*COEFFICIENT*inv;}
            float recip=num>ZERO?1.0f/num:INF,fac=factor(shapes[x][y][z][k],x,y,z,k),value;if(k==OUT){float recipO=numO>ZERO?1.0f/numO:INF;value=fac*((1.0f-oldOut)*recip+oldOut*recipO);out[x][y][z]=(fac*oldOut)*recipO;}else value=fac*recip;mem[x][y][z][k]=value;den+=value;objAt[k]=objectiveValue;}
            denMap[x][y][z]=den;if(den>0.0f)for(int k=0;k<K;k++){mem[x][y][z][k]/=den;if(k==OUT)out[x][y][z]/=den;obj+=objAt[k]*mem[x][y][z][k]*mem[x][y][z][k];}}
        objective[0]=obj;lastOutlier=out;lastDen=denMap;return mem;}
    private static int[] bits(float[][][][] a){int[] r=new int[NX*NY*NZ*a[0][0][0].length];int p=0;for(float[][][] q:a)for(float[][] w:q)for(float[] v:w)for(float z:v)r[p++]=Float.floatToIntBits(z);return r;}
    private static int[] bits(float[][][] a){int[] r=new int[NX*NY*NZ];int p=0;for(float[][] q:a)for(float[] w:q)for(float v:w)r[p++]=Float.floatToIntBits(v);return r;}
    private static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    public static void main(String[] args)throws Exception{images=readFloat4(args[0],NC);field=readFloat4(args[1],NC);bias=readFloat4(args[2],NC);initial=readFloat4(args[3],K);shapes=readFloat4(args[4],K);labels=readShort4(args[5]);initialOutlier=readFloat3(args[6]);initTables();StringBuilder s=new StringBuilder("{\"toads_source_sha256\":\"").append(TOADS_SHA).append("\",\"powertable_source_sha256\":\"").append(POWER_SHA).append("\",\"class_types\":[\"OPTIMIZED_GM\",\"OUTMASK\",\"MASK\",\"OBJECT\",\"OPTIMIZED_CSF\",\"OPTIMIZED_WM\"],\"cases\":[");boolean first=true;
        for(int f=0;f<2;f++)for(int b=0;b<2;b++){float[][][] auxEnergy=new float[NX][NY][NZ];float[][][][] energy=dataEnergy(f==1,b==1,auxEnergy);for(int sm=0;sm<2;sm++)for(int pr=0;pr<2;pr++){if(!first)s.append(',');first=false;float[] objective={0};float[][][][] mem=run(f==1,b==1,sm==1,pr==1,energy,auxEnergy,objective);s.append("{\"field\":").append(f==1).append(",\"psf\":").append(b==1).append(",\"smoothing\":").append(sm==1).append(",\"prior\":").append(pr==1).append(",\"class_energy_bits\":");append(s,bits(energy));s.append(",\"outlier_energy_bits\":");append(s,bits(auxEnergy));s.append(",\"membership_bits\":");append(s,bits(mem));s.append(",\"outlier_membership_bits\":");append(s,bits(lastOutlier));s.append(",\"denominator_bits\":");append(s,bits(lastDen));s.append(",\"objective_bits\":").append(Float.floatToIntBits(objective[0])).append('}');}}
        s.append("]}");System.out.println(s);}
}
