import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships excerpt in the default wrapper's class order. */
public final class ToadsDefaultAtlasMembershipReferenceHarness {
    private static final String TOADS_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWER_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int NX=9,NY=8,NZ=7,NC=2,K=11,MASK=0;
    private static final String[] TYPES={"MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_WM","OPTIMIZED_WM","OPTIMIZED_WM"};
    private static final float[][] CENTROID={{.28f,.31f,.34f,.37f,.40f,.43f,.46f,.49f,.52f,.55f,.58f},{.24f,.27f,.30f,.33f,.36f,.39f,.42f,.45f,.48f,.51f,.54f}};
    private static final float[][] VARIANCE={{.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},{1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}};
    private static final float[] MODALITY={.73f,.46f},IMAGE_VARIANCE={.61f,.34f};
    private static final float[] EXPCSF={.78f,.91f},EXPGM={1.16f,1.29f},EXPWM={1.55f,1.73f};
    private static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] SIM={
        {0,.11f,.14f,.18f,.13f,.16f,.12f,.15f,.09f,.10f,.08f},
        {.12f,0,.13f,.17f,.11f,.16f,.14f,.10f,.09f,.12f,.08f},
        {.15f,.12f,0,.13f,.18f,.12f,.16f,.11f,.10f,.07f,.09f},
        {.14f,.16f,.12f,0,.13f,.18f,.11f,.14f,.08f,.10f,.12f},
        {.11f,.13f,.17f,.12f,0,.11f,.19f,.14f,.10f,.09f,.07f},
        {.16f,.12f,.14f,.17f,.11f,0,.12f,.18f,.09f,.08f,.10f},
        {.13f,.15f,.11f,.14f,.18f,.12f,0,.13f,.08f,.09f,.11f},
        {.12f,.10f,.14f,.16f,.13f,.18f,.11f,0,.07f,.09f,.08f},
        {.10f,.09f,.12f,.08f,.11f,.09f,.10f,.07f,0,.13f,.16f},
        {.11f,.12f,.08f,.10f,.09f,.07f,.11f,.08f,.14f,0,.12f},
        {.08f,.10f,.09f,.12f,.07f,.11f,.10f,.09f,.17f,.13f,0}
    };
    private static final float OUTLIER=2.0f/K,SMOOTH=.043f,FACTOR=.64f,COEFFICIENT=.16f,ZERO=1e-10f,INF=1e10f;
    private static float[][][][] image,field,bias,initial,shapes;
    private static short[][][][] labels;
    private static int[][][] ribbon;
    private static PowerTable[][] powers;

    private static byte[] unzip(String path)throws Exception{ByteArrayOutputStream out=new ByteArrayOutputStream();try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] a=new byte[8192];int n;while((n=in.read(a))>=0)out.write(a,0,n);}return out.toByteArray();}
    private static float[][][][] readFloat4(String path,int nclass)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=nclass)throw new IllegalArgumentException("unexpected float NIfTI");int o=(int)b.getFloat(108);float[][][][] a=new float[NX][NY][NZ][nclass];for(int c=0;c<nclass;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getFloat(o+4*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static short[][][][] readShort4(String path)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=4)throw new IllegalArgumentException("unexpected int16 NIfTI");int o=(int)b.getFloat(108);short[][][][] a=new short[NX][NY][NZ][4];for(int c=0;c<4;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getShort(o+2*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static int[][][] readRibbon(String path)throws Exception{byte[] raw=unzip(path);ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(b.getInt(0)!=348||b.getShort(70)!=2||b.getShort(72)!=8||b.getShort(40)!=3||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ)throw new IllegalArgumentException("unexpected uint8 NIfTI");int o=(int)b.getFloat(108);int[][][] a=new int[NX][NY][NZ];for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z]=raw[o+(z*NY+y)*NX+x]&255;return a;}
    private static int cat(int k){if(TYPES[k].equals("OPTIMIZED_CSF"))return 0;if(TYPES[k].equals("OPTIMIZED_GM"))return 1;return 2;}
    private static float exponent(int category,int c){return category==0?EXPCSF[c]:category==1?EXPGM[c]:EXPWM[c];}
    private static void initTables(){powers=new PowerTable[3][NC];for(int t=0;t<3;t++)for(int c=0;c<NC;c++)powers[t][c]=new PowerTable(0.0f,4.0f,.000001f,exponent(t,c));}
    private static boolean direct(int x,int y,int z,int k){return labels[x][y][z][0]==k||labels[x][y][z][1]==k;}
    private static boolean neighbor(int x,int y,int z,int k){return labels[x][y][z][2]==k||labels[x][y][z][3]==k;}
    private static float shape(int x,int y,int z,int k,boolean prior){return prior?shapes[x][y][z][k]:0.0f;}
    private static float relation(int x,int y,int z,int k,boolean prior){float s=shape(x,y,z,k,prior);if(direct(x,y,z,k))return 1.0f;if(neighbor(x,y,z,k))return Math.max(.25f,Math.min(2.0f*FACTOR*s,1.0f));return Math.min(1.0f,2.0f*FACTOR*s);}
    private static float inversePrior(int x,int y,int z,int k){float v=0.0f;for(int m=0;m<K;m++)if(m!=k)v+=ATLAS_PRIOR[m]*SIM[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];return v;}
    private static int[] flatten(float[][][][] a){int[] r=new int[NX*NY*NZ*a[0][0][0].length];int p=0;for(float[][][] v:a)for(float[][] w:v)for(float[] q:w)for(float z:q)r[p++]=Float.floatToIntBits(z);return r;}
    private static int[] flatten(float[][][] a){int[] r=new int[NX*NY*NZ];int p=0;for(float[][] v:a)for(float[] w:v)for(float q:w)r[p++]=Float.floatToIntBits(q);return r;}
    private static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    private static float[][][][] makeData(boolean useField,boolean usePsf,boolean prior){float[][][][] out=new float[NX][NY][NZ][K];for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<K;k++){
        if(prior&&!direct(x,y,z,k)&&!neighbor(x,y,z,k)&&shapes[x][y][z][k]==0.0f)continue;
        if(!prior&&!direct(x,y,z,k)&&!neighbor(x,y,z,k))continue;
        float num=0.0f;if(TYPES[k].equals("MASK")){for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num+=w*(image[x][y][z][c]==0.0f?0.0f:1.0f);}}
        else{int t=cat(k);for(int c=0;c<NC;c++){float value=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];if(usePsf)value-=bias[x][y][z][c];float d=value-CENTROID[c][k];float dist=d*d,w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num=(float)(num+(double)w*powers[t][c].lookup((double)dist,(double)exponent(t,c)));}}
        out[x][y][z][k]=num;
    }return out;}
    private static float[][][][] run(boolean useField,boolean usePsf,boolean useSmooth,boolean prior,float[][][][] energy,float[][][] denMap,float[] objective){float[][][][] mem=new float[NX][NY][NZ][K];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++)mem[x][y][z][k]=initial[x][y][z][k];float obj=0.0f;int[][] ds={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){float den=0.0f;float[] objAt=new float[K];for(int k=0;k<K;k++){boolean d=direct(x,y,z,k),n=neighbor(x,y,z,k);float sh=shape(x,y,z,k,prior);if(!d&&!n&&sh==0.0f){mem[x][y][z][k]=0.0f;continue;}float num=energy[x][y][z][k],e=num;
            if(useSmooth){float ngb=0.0f;for(int[] v:ds)for(int m=0;m<K;m++)if(m!=k){float q=mem[x+v[0]][y+v[1]][z+v[2]][m];ngb+=q*q;}num+=SMOOTH_PRIOR[k]*SMOOTH*ngb;e+=.5f*SMOOTH_PRIOR[k]*SMOOTH*ngb;}
            if(prior){float inv=inversePrior(x,y,z,k);num+=COEFFICIENT*inv;e+=.5f*COEFFICIENT*inv;}
            float fac=relation(x,y,z,k,prior),recip=num>ZERO?1.0f/num:INF;mem[x][y][z][k]=fac*recip;den+=mem[x][y][z][k];objAt[k]=e;}
            denMap[x][y][z]=den;if(den>0)for(int k=0;k<K;k++){mem[x][y][z][k]/=den;obj+=objAt[k]*mem[x][y][z][k]*mem[x][y][z][k];}}
        objective[0]=obj;return mem;}
    public static void main(String[] args)throws Exception{image=readFloat4(args[0],NC);field=readFloat4(args[1],NC);bias=readFloat4(args[2],NC);initial=readFloat4(args[3],K);shapes=readFloat4(args[4],K);labels=readShort4(args[5]);ribbon=readRibbon(args[6]);initTables();StringBuilder s=new StringBuilder("{\"toads_source_sha256\":\"").append(TOADS_SHA).append("\",\"powertable_source_sha256\":\"").append(POWER_SHA).append("\",\"shape\":[9,8,7,11],\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],\"cases\":[");boolean first=true;
        for(int f=0;f<2;f++)for(int p=0;p<2;p++)for(int sm=0;sm<2;sm++)for(int pr=0;pr<2;pr++){if(!first)s.append(',');first=false;float[][][] den=new float[NX][NY][NZ];float[] obj={0};float[][][][] energy=makeData(f==1,p==1,pr==1);float[][][][] mem=run(f==1,p==1,sm==1,pr==1,energy,den,obj);s.append("{\"field\":").append(f==1).append(",\"psf\":").append(p==1).append(",\"smoothing\":").append(sm==1).append(",\"prior\":").append(pr==1).append(",\"class_energy_bits\":");append(s,flatten(energy));s.append(",\"membership_bits\":");append(s,flatten(mem));s.append(",\"denominator_bits\":");append(s,flatten(den));s.append(",\"objective_bits\":").append(Float.floatToIntBits(obj[0])).append('}');}
        s.append("]}");System.out.println(s);}
}
