import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships excerpt for default atlas classes plus OUTLIER. */
public final class ToadsDefaultOutlierMembershipReferenceHarness {
    private static final String TOADS_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWER_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int NX=9,NY=8,NZ=7,NC=2,BASE=11,K=12;
    private static final String[] BASE_TYPES={"MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_WM","OPTIMIZED_WM","OPTIMIZED_WM"};
    private static final float[][] CENTROID={{.28f,.31f,.34f,.37f,.40f,.43f,.46f,.49f,.52f,.55f,.58f},{.24f,.27f,.30f,.33f,.36f,.39f,.42f,.45f,.48f,.51f,.54f}};
    private static final float[][] VARIANCE={{.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},{1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}};
    private static final float[] MODALITY={.73f,.46f},IMAGE_VARIANCE={.61f,.34f};
    private static final float[] EXPCSF={.78f,.91f},EXPGM={1.16f,1.29f},EXPWM={1.55f,1.73f};
    private static final float[] SMOOTH_BASE={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_BASE={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] SIM_BASE={
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
    private static final float OUTLIER_VALUE=1.3f,OUT_CENTROID_0=.63f,OUT_CENTROID_1=.57f;
    private static final float OUT_VARIANCE_0=1.07f,OUT_VARIANCE_1=.94f,OUT_SMOOTH=.97f,OUT_ATLAS=.13f;
    private static final float SMOOTH=.043f,FACTOR=.64f,COEFFICIENT=.16f,ZERO=1e-10f,INF=1e10f;
    private static float[][][][] images,field,bias,baseInitial,baseShapes;
    private static float[][][] initialOutlier;
    private static short[][][][] baseLabels;
    private static int[][][] ribbon;
    private static PowerTable[][] powers;

    private static byte[] gunzip(String path)throws Exception{ByteArrayOutputStream out=new ByteArrayOutputStream();try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] block=new byte[8192];int n;while((n=in.read(block))>=0)out.write(block,0,n);}return out.toByteArray();}
    private static ByteBuffer nifti(String path)throws Exception{return ByteBuffer.wrap(gunzip(path)).order(ByteOrder.LITTLE_ENDIAN);}
    private static float[][][][] readFloat4(String path,int last)throws Exception{ByteBuffer b=nifti(path);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=last)throw new IllegalArgumentException("unexpected float NIfTI: "+path);int o=(int)b.getFloat(108);float[][][][] a=new float[NX][NY][NZ][last];for(int c=0;c<last;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getFloat(o+4*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static float[][][] readFloat3(String path)throws Exception{ByteBuffer b=nifti(path);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=3||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ)throw new IllegalArgumentException("unexpected float NIfTI: "+path);int o=(int)b.getFloat(108);float[][][] a=new float[NX][NY][NZ];for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z]=b.getFloat(o+4*((z*NY+y)*NX+x));return a;}
    private static short[][][][] readShort4(String path)throws Exception{ByteBuffer b=nifti(path);if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ||b.getShort(48)!=4)throw new IllegalArgumentException("unexpected labels NIfTI");int o=(int)b.getFloat(108);short[][][][] a=new short[NX][NY][NZ][4];for(int c=0;c<4;c++)for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z][c]=b.getShort(o+2*(((c*NZ+z)*NY+y)*NX+x));return a;}
    private static int[][][] readRibbon(String path)throws Exception{ByteBuffer b=nifti(path);if(b.getInt(0)!=348||b.getShort(70)!=2||b.getShort(72)!=8||b.getShort(40)!=3||b.getShort(42)!=NX||b.getShort(44)!=NY||b.getShort(46)!=NZ)throw new IllegalArgumentException("unexpected ribbon NIfTI");int o=(int)b.getFloat(108);byte[] raw=b.array();int[][][] a=new int[NX][NY][NZ];for(int z=0;z<NZ;z++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++)a[x][y][z]=raw[o+(z*NY+y)*NX+x]&255;return a;}
    private static int baseIndex(int k,int outIndex){return k<outIndex?k:k==outIndex?-1:k-1;}
    private static String type(int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?"OUTLIER":BASE_TYPES[b];}
    private static float centroid(int c,int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?(c==0?OUT_CENTROID_0:OUT_CENTROID_1):CENTROID[c][b];}
    private static float variance(int c,int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?(c==0?OUT_VARIANCE_0:OUT_VARIANCE_1):VARIANCE[c][b];}
    private static float smoothPrior(int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?OUT_SMOOTH:SMOOTH_BASE[b];}
    private static float atlasPrior(int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?OUT_ATLAS:ATLAS_BASE[b];}
    private static float similarity(int k,int m,int outIndex){int a=baseIndex(k,outIndex),b=baseIndex(m,outIndex);if(a>=0&&b>=0)return SIM_BASE[a][b];if(k==m)return 0.0f;return .05f+((2*k+3*m)%7)*.01f;}
    private static float initial(int x,int y,int z,int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?.16f+((x+2*y+z)%7)*.014f:baseInitial[x][y][z][b];}
    private static float shape(int x,int y,int z,int k,int outIndex){int b=baseIndex(k,outIndex);return b<0?.29f+.02f*ribbon[x][y][z]+((x+2*y+z)%5)*.004f:baseShapes[x][y][z][b];}
    private static int label(int x,int y,int z,int slot,int outIndex){if(slot==0&&(x+2*y+z)%3==0)return outIndex;if(slot==2&&(x+2*y+z)%3==1)return outIndex;int value=baseLabels[x][y][z][slot];return value<0?value:(value>=outIndex?value+1:value);}
    private static boolean direct(int x,int y,int z,int k,int outIndex){return label(x,y,z,0,outIndex)==k||label(x,y,z,1,outIndex)==k;}
    private static boolean neighbor(int x,int y,int z,int k,int outIndex){return label(x,y,z,2,outIndex)==k||label(x,y,z,3,outIndex)==k;}
    private static float relationShape(int x,int y,int z,int k,int outIndex,boolean prior){return prior?shape(x,y,z,k,outIndex):0.0f;}
    private static float relation(int x,int y,int z,int k,int outIndex,boolean prior){float s=relationShape(x,y,z,k,outIndex,prior);if(direct(x,y,z,k,outIndex))return 1.0f;if(neighbor(x,y,z,k,outIndex))return Math.max(.25f,Math.min(2.0f*FACTOR*s,1.0f));return Math.min(1.0f,2.0f*FACTOR*s);}
    private static int cat(String kind){return kind.equals("OPTIMIZED_CSF")?0:kind.equals("OPTIMIZED_GM")?1:2;}
    private static float exponent(int category,int c){return category==0?EXPCSF[c]:category==1?EXPGM[c]:EXPWM[c];}
    private static void initTables(){powers=new PowerTable[3][NC];for(int t=0;t<3;t++)for(int c=0;c<NC;c++)powers[t][c]=new PowerTable(0.0f,4.0f,.000001f,exponent(t,c));}
    private static float dataTerm(int x,int y,int z,int k,int outIndex,boolean useField,boolean usePsf){String kind=type(k,outIndex);float value=0.0f;if(kind.equals("OUTLIER")){for(int c=0;c<NC;c++){float corrected=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];if(usePsf)corrected-=bias[x][y][z][c];float residual=corrected-centroid(c,k,outIndex);float w=MODALITY[c]*variance(c,k,outIndex)*IMAGE_VARIANCE[c];value+=w*residual*residual;}return value;}
        if(kind.equals("MASK")){for(int c=0;c<NC;c++){float w=MODALITY[c]*variance(c,k,outIndex)*IMAGE_VARIANCE[c];value+=w*(images[x][y][z][c]==0.0f?0.0f:1.0f);}return value;}
        int category=cat(kind);for(int c=0;c<NC;c++){float corrected=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];if(usePsf)corrected-=bias[x][y][z][c];float residual=corrected-centroid(c,k,outIndex);float distance=residual*residual;float w=MODALITY[c]*variance(c,k,outIndex)*IMAGE_VARIANCE[c];value=(float)(value+(double)w*powers[category][c].lookup((double)distance,(double)exponent(category,c)));}return value;}
    private static float auxTerm(int k,int outIndex){float value=0.0f;for(int c=0;c<NC;c++){float w=MODALITY[c]*variance(c,k,outIndex);w=w*IMAGE_VARIANCE[c];value+=w*OUTLIER_VALUE;}return value;}
    private static float[][][][] insertBase(float[][][][] base,int outIndex){float[][][][] result=new float[NX][NY][NZ][K];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){int b=baseIndex(k,outIndex);result[x][y][z][k]=b<0?0.0f:base[x][y][z][b];}return result;}
    private static float[][][][] run(int outIndex,boolean useField,boolean usePsf,boolean useSmooth,boolean usePrior,float[][][][] energy,float[][][] aux,float[][][] denominator,float[] objective){String[] types=new String[K];float[] sm=new float[K],at=new float[K];float[][] sim=new float[K][K];float[][][][] memberships=insertBase(baseInitial,outIndex),shapes=insertBase(baseShapes,outIndex);float[][][] outMap=new float[NX][NY][NZ];
        for(int k=0;k<K;k++){types[k]=type(k,outIndex);sm[k]=smoothPrior(k,outIndex);at[k]=atlasPrior(k,outIndex);for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){memberships[x][y][z][k]=initial(x,y,z,k,outIndex);shapes[x][y][z][k]=shape(x,y,z,k,outIndex);}for(int m=0;m<K;m++)sim[k][m]=similarity(k,m,outIndex);}
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)outMap[x][y][z]=initialOutlier[x][y][z];
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<K;k++){
            if(usePrior&&!direct(x,y,z,k,outIndex)&&!neighbor(x,y,z,k,outIndex)&&shape(x,y,z,k,outIndex)==0.0f)continue;
            if(!usePrior&&!direct(x,y,z,k,outIndex)&&!neighbor(x,y,z,k,outIndex))continue;
            float num=dataTerm(x,y,z,k,outIndex,useField,usePsf);energy[x][y][z][k]=num;
            if(k==outIndex){float numO=auxTerm(k,outIndex);aux[x][y][z]=numO;float out=outMap[x][y][z];energy[x][y][z][k]=(1.0f-out)*num+out*numO;}
        }
        float objectiveTotal=0.0f;int[][] ds={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){float den=0.0f;float[] objectiveAt=new float[K];for(int k=0;k<K;k++){
            if(!direct(x,y,z,k,outIndex)&&!neighbor(x,y,z,k,outIndex)&&relationShape(x,y,z,k,outIndex,usePrior)==0.0f){memberships[x][y][z][k]=0.0f;continue;}
            float num=energy[x][y][z][k],numO=aux[x][y][z],oldOut=outMap[x][y][z],obj=num;
            if(useSmooth){float ngb=0.0f;for(int[] d:ds)for(int m=0;m<K;m++)if(m!=k){float q=memberships[x+d[0]][y+d[1]][z+d[2]][m];ngb+=q*q;}float add=sm[k]*SMOOTH*ngb;num+=add;obj+=(.5f*sm[k])*SMOOTH*ngb;}
            if(usePrior){float inv=0.0f;for(int m=0;m<K;m++)if(m!=k){float term=at[m]*sim[k][m];term=term*shapes[x][y][z][m];term=term*shapes[x][y][z][m];inv+=term;}float add=COEFFICIENT*inv;num+=add;if(k==outIndex)numO+=add;obj+=(.5f*COEFFICIENT)*inv;}
            float reciprocal=num>ZERO?1.0f/num:INF,value;float rel=relation(x,y,z,k,outIndex,usePrior);
            if(k==outIndex){float reciprocalO=numO>ZERO?1.0f/numO:INF;value=rel*((1.0f-oldOut)*reciprocal+oldOut*reciprocalO);outMap[x][y][z]=(rel*oldOut)*reciprocalO;}else value=rel*reciprocal;
            memberships[x][y][z][k]=value;den+=value;objectiveAt[k]=obj;
        }denominator[x][y][z]=den;if(den>0.0f)for(int k=0;k<K;k++){memberships[x][y][z][k]=memberships[x][y][z][k]/den;if(k==outIndex)outMap[x][y][z]=outMap[x][y][z]/den;objectiveTotal+=objectiveAt[k]*memberships[x][y][z][k]*memberships[x][y][z][k];}}
        objective[0]=objectiveTotal;lastOutMap=outMap;return memberships;}
    private static int[] bits(float[][][][] a){int[] r=new int[NX*NY*NZ*a[0][0][0].length];int p=0;for(float[][][] q:a)for(float[][] v:q)for(float[] w:v)for(float f:w)r[p++]=Float.floatToIntBits(f);return r;}
    private static int[] bits(float[][][] a){int[] r=new int[NX*NY*NZ];int p=0;for(float[][] v:a)for(float[] w:v)for(float f:w)r[p++]=Float.floatToIntBits(f);return r;}
    private static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    private static List<int[]> cases(){List<int[]> result=new ArrayList<>();for(int idx:new int[]{0,6,11})result.add(new int[]{idx,1,1,1,1});for(int f=0;f<2;f++)for(int p=0;p<2;p++)for(int sm=0;sm<2;sm++)for(int pr=0;pr<2;pr++)result.add(new int[]{6,f,p,sm,pr});return result;}
    public static void main(String[] args)throws Exception{images=readFloat4(args[0],NC);field=readFloat4(args[1],NC);bias=readFloat4(args[2],NC);baseInitial=readFloat4(args[3],BASE);baseShapes=readFloat4(args[4],BASE);baseLabels=readShort4(args[5]);ribbon=readRibbon(args[6]);initialOutlier=readFloat3(args[7]);initTables();StringBuilder s=new StringBuilder("{\"toads_source_sha256\":\"").append(TOADS_SHA).append("\",\"powertable_source_sha256\":\"").append(POWER_SHA).append("\",\"shape\":[9,8,7,12],\"default_class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],\"cases\":[");boolean first=true;for(int[] c:cases()){if(!first)s.append(',');first=false;int outIndex=c[0];float[][][][] energy=new float[NX][NY][NZ][K];float[][][] aux=new float[NX][NY][NZ],den=new float[NX][NY][NZ];float[] obj={0};float[][][][] mem=run(outIndex,c[1]==1,c[2]==1,c[3]==1,c[4]==1,energy,aux,den,obj);s.append("{\"outlier_index\":").append(outIndex).append(",\"field\":").append(c[1]==1).append(",\"psf\":").append(c[2]==1).append(",\"smoothing\":").append(c[3]==1).append(",\"prior\":").append(c[4]==1).append(",\"types\":[");for(int k=0;k<K;k++){if(k>0)s.append(',');s.append('"').append(type(k,outIndex)).append('"');}s.append("],\"class_energy_bits\":");append(s,bits(energy));s.append(",\"outlier_energy_bits\":");append(s,bits(aux));s.append(",\"membership_bits\":");append(s,bits(mem));s.append(",\"outlier_map_bits\":");append(s,bits(lastOutMap));s.append(",\"denominator_bits\":");append(s,bits(den));s.append(",\"objective_bits\":").append(Float.floatToIntBits(obj[0])).append('}');}s.append("]}");System.out.println(s);}
    private static float[][][] lastOutMap;
}
