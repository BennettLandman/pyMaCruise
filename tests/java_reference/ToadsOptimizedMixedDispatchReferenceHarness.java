import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected computeExactMemberships excerpt for interleaved OBJECT/MASK/optimized classes. */
public final class ToadsOptimizedMixedDispatchReferenceHarness {
    private static final String TOADS_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWERTABLE_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int NX=15, NY=13, NZ=11, CHANNELS=3, CLASSES=5;
    private static final String[] TYPES={"OBJECT","OPTIMIZED_CSF","MASK","OPTIMIZED_GM","OPTIMIZED_WM"};
    private static final float[][] CENTROID={
        {0.30f,0.35f,0.40f,0.45f,0.50f},
        {0.26f,0.31f,0.36f,0.41f,0.46f},
        {0.22f,0.27f,0.32f,0.37f,0.42f}
    };
    private static final float[][] VARIANCE={
        {0.8f,0.9f,1.0f,1.1f,1.2f},
        {1.1f,0.7f,0.9f,1.3f,0.8f},
        {0.6f,1.2f,0.75f,0.95f,1.15f}
    };
    private static final float[] MODALITY={0.7f,0.5f,0.3f};
    private static final float[] IMAGE_VARIANCE={0.55f,0.30f,0.15f};
    // Indexed by optimized tissue family: CSF, GM, WM.
    private static final float[] EXPONENT={0.8f,1.2f,1.6f};
    private static final float[] SMOOTH_PRIOR={0.5f,0.7f,0.9f,1.1f,1.3f};
    private static final float[] ATLAS_PRIOR={0.12f,0.18f,0.22f,0.25f,0.23f};
    private static final float[][] SIMILARITY={
        {0,.12f,.18f,.21f,.16f}, {.15f,0,.13f,.19f,.17f},
        {.18f,.14f,0,.11f,.20f}, {.16f,.18f,.13f,0,.12f},
        {.19f,.15f,.17f,.14f,0}
    };
    private static final float SMOOTH=0.027f, FACTOR=0.64f, COEFFICIENT=0.11f;
    private static final float ZERO=1e-10f, INF=1e10f;
    private static float[][][][] images, field, bias, memberships, shapes, dataTerms;
    private static int[][][] ribbon;
    private static int[][][][] closest;
    private static final PowerTable[] POWERS={
        new PowerTable(0.0f,4.0f,0.000001f,0.8f),
        new PowerTable(0.0f,4.0f,0.000001f,1.2f),
        new PowerTable(0.0f,4.0f,0.000001f,1.6f)
    };

    private static byte[] gzip(String path) throws IOException {
        ByteArrayOutputStream out=new ByteArrayOutputStream();
        try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){
            byte[] buffer=new byte[8192]; int n;
            while((n=in.read(buffer))>=0) out.write(buffer,0,n);
        }
        return out.toByteArray();
    }
    private static float[][][][] read4d(String path) throws IOException {
        byte[] raw=gzip(path); ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348 || b.getShort(70)!=16 || b.getShort(72)!=32 || b.getShort(40)!=4 ||
                b.getShort(42)!=NX || b.getShort(44)!=NY || b.getShort(46)!=NZ || b.getShort(48)!=CHANNELS)
            throw new IOException("unexpected float32 fixture dimensions/type");
        int offset=(int)b.getFloat(108); float[][][][] a=new float[NX][NY][NZ][CHANNELS];
        for(int c=0;c<CHANNELS;c++) for(int z=0;z<NZ;z++) for(int y=0;y<NY;y++) for(int x=0;x<NX;x++){
            int i=(((c*NZ+z)*NY+y)*NX+x); a[x][y][z][c]=b.getFloat(offset+4*i);
        }
        return a;
    }
    private static int[][][] readRibbon(String path) throws IOException {
        byte[] raw=gzip(path); ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        if(b.getInt(0)!=348 || b.getShort(70)!=2 || b.getShort(72)!=8 || b.getShort(40)!=3 ||
                b.getShort(42)!=NX || b.getShort(44)!=NY || b.getShort(46)!=NZ)
            throw new IOException("unexpected uint8 ribbon fixture");
        int offset=(int)b.getFloat(108); int[][][] a=new int[NX][NY][NZ];
        for(int z=0;z<NZ;z++) for(int y=0;y<NY;y++) for(int x=0;x<NX;x++){
            int i=(z*NY+y)*NX+x; a[x][y][z]=raw[offset+i]&255;
        }
        return a;
    }
    private static void initialize(){
        memberships=new float[NX][NY][NZ][CLASSES]; shapes=new float[NX][NY][NZ][CLASSES];
        closest=new int[NX][NY][NZ][4];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++){
            for(int k=0;k<CLASSES;k++){
                int im=k==0?(x+y+z)%5:(x+2*y+z+3*k)%7;
                memberships[x][y][z][k]=0.12f+im*0.017f;
                int n=(x+2*y+3*z+k)%5;
                shapes[x][y][z][k]=0.2f+0.5f*ribbon[x][y][z]+n*0.01f;
            }
            java.util.Arrays.fill(closest[x][y][z],-1);
            closest[x][y][z][0]=(x+y+z)%CLASSES;
            if((x+2*y+z)%7==0) closest[x][y][z][1]=(closest[x][y][z][0]+2)%CLASSES;
            if((x+y+2*z)%3==0) closest[x][y][z][2]=(closest[x][y][z][0]+1)%CLASSES;
            if((x+z)%4==0) closest[x][y][z][3]=(closest[x][y][z][0]+3)%CLASSES;
        }
    }
    private static int[] update(boolean useField,boolean usePsf){
        initialize(); dataTerms=new float[NX][NY][NZ][CLASSES]; float total=0.0f;
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++){
            float den=0; float[] objective=new float[CLASSES];
            for(int k=0;k<CLASSES;k++){
                float num=0;
                for(int c=0;c<CHANNELS;c++){
                    float corrected=useField?field[x][y][z][c]*images[x][y][z][c]:images[x][y][z][c];
                    if(usePsf) corrected-=bias[x][y][z][c];
                    float weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];
                    if(TYPES[k].equals("OBJECT")){
                        float residual=corrected-CENTROID[c][k];
                        num+=weight*residual*residual;
                    }else if(TYPES[k].equals("MASK")){
                        num+=weight*(images[x][y][z][c]==0.0f?0.0f:1.0f);
                    }else{
                        float residual=corrected-CENTROID[c][k]; float distance=residual*residual;
                        int t=TYPES[k].equals("OPTIMIZED_CSF")?0:TYPES[k].equals("OPTIMIZED_GM")?1:2;
                        double powered=POWERS[t].lookup(distance,EXPONENT[t]);
                        num+=(double)weight*powered;
                    }
                }
                dataTerms[x][y][z][k]=num;
                float energy=num;
                if(SMOOTH>0){
                    float neighbor=0;
                    for(int[] d:offsets) for(int m=0;m<CLASSES;m++) if(m!=k){float v=memberships[x+d[0]][y+d[1]][z+d[2]][m];neighbor+=v*v;}
                    num+=SMOOTH_PRIOR[k]*SMOOTH*neighbor;
                    energy+=0.5f*SMOOTH_PRIOR[k]*SMOOTH*neighbor;
                }
                float invprior=0;
                for(int m=0;m<CLASSES;m++) if(m!=k){
                    float term=ATLAS_PRIOR[m]*SIMILARITY[k][m];term*=shapes[x][y][z][m];term*=shapes[x][y][z][m];invprior+=term;
                }
                num+=COEFFICIENT*invprior; energy+=0.5f*COEFFICIENT*invprior;
                boolean direct=closest[x][y][z][0]==k||closest[x][y][z][1]==k;
                boolean neighborClass=closest[x][y][z][2]==k||closest[x][y][z][3]==k;
                float factor;
                if(direct) factor=1.0f;
                else {float scaled=2.0f*FACTOR*shapes[x][y][z][k];factor=neighborClass?Math.max(0.25f,Math.min(scaled,1.0f)):Math.min(1.0f,scaled);}
                float reciprocal=num>ZERO?1.0f/num:INF;
                memberships[x][y][z][k]=factor*reciprocal; den+=memberships[x][y][z][k]; objective[k]=energy;
            }
            if(den>0) for(int k=0;k<CLASSES;k++){
                memberships[x][y][z][k]/=den;
                total+=objective[k]*memberships[x][y][z][k]*memberships[x][y][z][k];
            }
        }
        int[] result=new int[NX*NY*NZ*CLASSES+1];int p=0;
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int k=0;k<CLASSES;k++) result[p++]=Float.floatToIntBits(memberships[x][y][z][k]);
        result[p]=Float.floatToIntBits(total); return result;
    }
    private static void printBits(StringBuilder out,int[] bits){out.append('[');for(int i=0;i<bits.length;i++){if(i>0)out.append(',');out.append(bits[i]);}out.append(']');}
    private static int[] bits(float[][][][] a){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(a[x][y][z][k]);return b;}
    public static void main(String[] args)throws Exception{
        images=read4d(args[0]);field=read4d(args[1]);bias=read4d(args[2]);ribbon=readRibbon(args[3]);
        boolean[][] modes={{false,false},{false,true},{true,false},{true,true}};
        System.out.print("{\"toads_source_sha256\":\""+TOADS_SHA+"\",\"powertable_source_sha256\":\""+POWERTABLE_SHA+"\",\"shape\":[15,13,11,5],\"class_types\":[\"OBJECT\",\"OPTIMIZED_CSF\",\"MASK\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\"],\"cases\":[");
        for(int i=0;i<modes.length;i++){
            if(i>0)System.out.print(',');int[] bits=update(modes[i][0],modes[i][1]);
            int[] membershipsOnly=java.util.Arrays.copyOf(bits,bits.length-1);
            StringBuilder b=new StringBuilder(),energyBits=new StringBuilder();printBits(b,membershipsOnly);printBits(energyBits,bits(dataTerms));
            System.out.print("{\"field\":"+modes[i][0]+",\"psf\":"+modes[i][1]+",\"updated_bits\":"+b+
                ",\"data_term_bits\":"+energyBits+",\"objective_bits\":"+bits[bits.length-1]+"}");
        }
        System.out.println("]}");
    }
}
