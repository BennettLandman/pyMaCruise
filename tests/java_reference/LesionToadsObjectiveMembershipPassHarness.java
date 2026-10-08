import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.Arrays;

/**
 * P39 objective membership excerpt: OBJECT + LESION + BLACKHOLE + OUTMASK,
 * including the selected zero-shape class skip and stale energy behavior.
 * Source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * The synthetic arrays use the same exact binary-fraction fixture generator as
 * tests/data/synthetic/lesion_objective_membership_pass/generate_fixture.py.
 */
public final class LesionToadsObjectiveMembershipPassHarness {
    private static final int NX=5, NY=4, NZ=5, K=4;
    private static final String[] TYPES={"OBJECT","LESION","BLACKHOLE","OUTMASK"};
    private static final boolean[] PRIOR={false,true,false,true};
    private static final float[] SMOOTH_PRIOR={0.375f,0.5f,0.625f,0.75f};
    private static final float[] ATLAS_PRIOR={0.125f,0.25f,0.375f,0.5f};
    private static final float[][] SIMILARITY={
        {0.0f,0.3125f,0.4375f,0.5625f},
        {0.25f,0.0f,0.5f,0.625f},
        {0.375f,0.3125f,0.0f,0.6875f},
        {0.1875f,0.4375f,0.5625f,0.0f}
    };
    private static final float[] LESION_SIM={0.6875f,0.5625f,0.4375f,0.3125f};
    private static final float[] BLACKHOLE_SIM={0.25f,0.375f,0.5f,0.625f};
    private static final float SMOOTHING=0.25f, PRIOR_COEFFICIENT=0.375f, ZERO=1.0e-10f, INF=1.0e10f;
    private static int bits(float f){return Float.floatToIntBits(f);}

    private static float[][][][] memberships(){
        float[][][][] a=new float[NX][NY][NZ][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int fold=(x*7+y*11+z*5+x*y)%29;
            for(int k=0;k<K;k++)a[x][y][z][k]=((fold+3*k)%17+1)/32.0f;
        } return a;
    }
    private static float[][][] scalar(int offset,int multiplier,int modulus,int base,int denominator){
        float[][][] a=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int fold=(x*7+y*11+z*5+x*y)%29;
            a[x][y][z]=((fold*multiplier)%modulus+base)/(float)denominator;
        } return a;
    }
    private static float[][][][] classVolume(int firstChannel){
        float[][][][] a=new float[NX][NY][NZ][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int fold=(x*7+y*11+z*5+x*y)%29;
            for(int k=0;k<K;k++){
                int value;
                if(firstChannel==7)value=(fold*3+5*k)%23+4;
                else value=(x*3+y*5+z*7+k*2)%12+4;
                a[x][y][z][k]=value/(float)(firstChannel==7?16:16);
            }
        } return a;
    }
    private static float[] classShape(float[][][][] all,int x,int y,int z){return all[x][y][z];}
    private static float shapeFactor(float v){float doubled=2.0f*v;return 1.0f<doubled?1.0f:doubled;}
    private static float invert(float value){return value>ZERO?1.0f/value:INF;}
    private static float[][][] copy3(float[][][] in){float[][][] o=new float[NX][NY][NZ];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)o[x][y][z]=in[x][y][z];return o;}
    private static float[][][][] copy4(float[][][][] in){float[][][][] o=new float[NX][NY][NZ][K];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)o[x][y][z]=in[x][y][z].clone();return o;}

    private static final class Result {
        float[][][][] mem,energy; float[][][] lesion,bh,outlier,den,lesionEnergy,bhEnergy,running; float objective;
    }

    private static Result run(){
        float[][][][] mem=memberships();
        float[][][] lesion=scalar(4,5,13,2,16), bh=scalar(5,7,11,2,16), outlier=scalar(6,2,7,1,16);
        float[][][][] data=classVolume(7), shapes=classVolume(14);
        shapes[2][1][2][0]=0.0f;
        shapes[2][1][3][1]=0.0f;
        shapes[2][2][1][2]=0.0f;
        shapes[2][2][2][3]=0.0f;
        float[][][] lesionData=scalar(11,9,19,3,16), bhData=scalar(12,4,17,2,16), outmaskData=scalar(13,6,15,2,16);
        float[][][][] energy=new float[NX][NY][NZ][K];
        float[][][] lesionEnergy=new float[NX][NY][NZ], bhEnergy=new float[NX][NY][NZ], denMap=new float[NX][NY][NZ];
        float[] energyState=new float[K]; float lesionEnergyState=0.0f, bhEnergyState=0.0f;
        float total=0.0f; float[][][] running=new float[NX][NY][NZ];
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0.0f;
            for(int k=0;k<K;k++){
                if(shapes[x][y][z][k]<=0.0f){
                    // The selected Java loop sets skipped class membership
                    // to zero after the shape[k] > 0 branch, while energy[k]
                    // is left at its last value for objective calculation.
                    mem[x][y][z][k]=0.0f;
                    energy[x][y][z][k]=energyState[k];
                    continue;
                }
                float num=data[x][y][z][k], numL=0.0f, numB=0.0f, numO=0.0f;
                if(k==1||k==2)numL=lesionData[x][y][z];
                if(k==2)numB=bhData[x][y][z];
                if(k==3)numO=outmaskData[x][y][z];
                float lesionE=(k==1||k==2)?numL:0.0f, bhE=(k==2)?numB:0.0f;
                float e=(k==3)?(1.0f-outlier[x][y][z])*num+outlier[x][y][z]*numO:num;
                float[] neighbors=new float[6]; int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                float ngb=0,ngbL=0,ngbB=0;
                for(int f=0;f<6;f++){
                    int xx=x+delta[f][0],yy=y+delta[f][1],zz=z+delta[f][2];
                    for(int m=0;m<K;m++)if(m!=k)ngb+=mem[xx][yy][zz][m]*mem[xx][yy][zz][m];
                    ngbL+=lesion[xx][yy][zz]*lesion[xx][yy][zz];
                    ngbB+=bh[xx][yy][zz]*bh[xx][yy][zz];
                }
                float ngbS=0;
                if(k==1||k==2){
                    int[][] self={{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
                    for(int[] d:self)ngbS+=mem[x+d[0]][y+d[1]][z+d[2]][k]*mem[x+d[0]][y+d[1]][z+d[2]][k];
                }
                float p=SMOOTH_PRIOR[k];
                num+=p*SMOOTHING*(ngb+ngbL+ngbB);
                e+=0.5f*p*SMOOTHING*(ngb+ngbL+ngbB);
                if(k==1||k==2){numL+=p*SMOOTHING*(ngb+ngbS+ngbB);}
                if(k==2){numB+=p*SMOOTHING*(ngb+ngbS+ngbL);}
                if(k==1||k==2)lesionE+=0.5f*p*SMOOTHING*(ngb+ngbS+ngbB);
                if(k==2)bhE+=0.5f*p*SMOOTHING*(ngb+ngbS+ngbL);
                if(k==3)numO+=p*SMOOTHING*(ngb+ngbL);
                float inv=0,invL=0,invB=0;
                for(int m=0;m<K;m++){
                    if(m!=k)inv+=ATLAS_PRIOR[m]*SIMILARITY[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];
                    if(k==1||k==2)invL+=ATLAS_PRIOR[m]*LESION_SIM[m]*shapes[x][y][z][m]*shapes[x][y][z][m];
                    if(k==2)invB+=ATLAS_PRIOR[m]*BLACKHOLE_SIM[m]*shapes[x][y][z][m]*shapes[x][y][z][m];
                }
                inv+=ATLAS_PRIOR[2]*LESION_SIM[k]*shapes[x][y][z][2]*shapes[x][y][z][2];
                inv+=ATLAS_PRIOR[2]*BLACKHOLE_SIM[k]*shapes[x][y][z][2]*shapes[x][y][z][2];
                num+=PRIOR_COEFFICIENT*inv;
                if(k==1||k==2)numL+=PRIOR_COEFFICIENT*invL;
                if(k==2)numB+=PRIOR_COEFFICIENT*invB;
                if(k==3)numO+=PRIOR_COEFFICIENT*inv;
                e+=0.5f*PRIOR_COEFFICIENT*inv;
                if(k==1||k==2)lesionE+=0.5f*PRIOR_COEFFICIENT*invL;
                if(k==2)bhE+=0.5f*PRIOR_COEFFICIENT*invB;
                energyState[k]=e;
                if(k==1||k==2)lesionEnergyState=lesionE;
                if(k==2)bhEnergyState=bhE;
                energy[x][y][z][k]=energyState[k];
                float classMem=invert(num), lesionMem=0, bhMem=0;
                if(k==1||k==2)lesionMem=invert(numL);
                if(k==2)bhMem=invert(numB);
                if(PRIOR[k]){
                    float factor=shapeFactor(shapes[x][y][z][k]); classMem*=factor;
                    if(k==1||k==2)lesionMem*=factor;
                    if(k==2)bhMem*=factor;
                }
                if(k==3){
                    float classOut=invert(num), outMem=invert(numO);
                    float factor=PRIOR[k]?shapeFactor(shapes[x][y][z][k]):1.0f;
                    mem[x][y][z][k]=factor*((1.0f-outlier[x][y][z])*classOut+outlier[x][y][z]*outMem);
                    outlier[x][y][z]=factor*outlier[x][y][z]*outMem;
                    classMem=mem[x][y][z][k];
                }else mem[x][y][z][k]=classMem;
                if(k==1||k==2)lesion[x][y][z]=lesionMem;
                if(k==2)bh[x][y][z]=bhMem;
                if(k==1)den+=lesion[x][y][z];
                else if(k==2)den+=lesion[x][y][z]+bh[x][y][z];
                den+=mem[x][y][z][k];
            }
            denMap[x][y][z]=den;
            lesionEnergy[x][y][z]=lesionEnergyState;
            bhEnergy[x][y][z]=bhEnergyState;
            if(den>0.0f){
                for(int k=0;k<K;k++){
                    mem[x][y][z][k]/=den;
                    if(k==1||k==2)lesion[x][y][z]/=den;
                    if(k==2)bh[x][y][z]/=den;
                    else if(k==3)outlier[x][y][z]/=den;
                }
                for(int k=0;k<K;k++)total+=energy[x][y][z][k]*mem[x][y][z][k]*mem[x][y][z][k];
            }
            total+=lesionEnergy[x][y][z]*lesion[x][y][z]*lesion[x][y][z]+bhEnergy[x][y][z]*bh[x][y][z]*bh[x][y][z];
            running[x][y][z]=total;
        }
        Result r=new Result();r.mem=mem;r.energy=energy;r.lesion=lesion;r.bh=bh;r.outlier=outlier;r.den=denMap;r.lesionEnergy=lesionEnergy;r.bhEnergy=bhEnergy;r.running=running;r.objective=total;return r;
    }

    private static String bits(float[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){if(!first)s.append(',');first=false;s.append(bits(a[x][y][z][k]));}return s.append(']').toString();}
    private static String bits(float[][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){if(!first)s.append(',');first=false;s.append(bits(a[x][y][z]));}return s.append(']').toString();}
    private static String inputDigest() throws Exception {
        MessageDigest md=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        float[][][][] mem=memberships(), data=classVolume(7), shapes=classVolume(14);
        shapes[2][1][2][0]=0.0f;
        shapes[2][1][3][1]=0.0f;
        shapes[2][2][1][2]=0.0f;
        shapes[2][2][2][3]=0.0f;
        float[][][] a=scalar(4,5,13,2,16), c=scalar(5,7,11,2,16), d=scalar(6,2,7,1,16), e=scalar(11,9,19,3,16), f=scalar(12,4,17,2,16), g=scalar(13,6,15,2,16);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            for(int k=0;k<K;k++){put(md,b,mem[x][y][z][k]);}
            put(md,b,a[x][y][z]);put(md,b,c[x][y][z]);put(md,b,d[x][y][z]);
            for(int k=0;k<K;k++)put(md,b,data[x][y][z][k]);
            put(md,b,e[x][y][z]);put(md,b,f[x][y][z]);put(md,b,g[x][y][z]);
            for(int k=0;k<K;k++)put(md,b,shapes[x][y][z][k]);
        }
        StringBuilder out=new StringBuilder();for(byte v:md.digest())out.append(String.format("%02x",v));return out.toString();
    }
    private static void put(MessageDigest md,ByteBuffer b,float v){b.clear();b.putInt(bits(v));md.update(b.array());}
    public static void main(String[] args)throws Exception{
        Result r=run();
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"input_sha256\":\""+inputDigest()+"\",\"shape\":["+NX+","+NY+","+NZ+"],\"class_types\":[\"OBJECT\",\"LESION\",\"BLACKHOLE\",\"OUTMASK\"],\"class_memberships_bits\":"+bits(r.mem)+",\"lesion_membership_bits\":"+bits(r.lesion)+",\"black_hole_membership_bits\":"+bits(r.bh)+",\"outlier_membership_bits\":"+bits(r.outlier)+",\"denominator_bits\":"+bits(r.den)+",\"class_energy_bits\":"+bits(r.energy)+",\"lesion_energy_bits\":"+bits(r.lesionEnergy)+",\"black_hole_energy_bits\":"+bits(r.bhEnergy)+",\"running_objective_bits\":"+bits(r.running)+",\"objective_bits\":"+bits(r.objective)+"}");
    }
}
