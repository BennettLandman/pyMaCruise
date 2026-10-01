import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.DataInputStream;
import java.io.FileInputStream;

/** Fixed-input source-order standard-13 membership pass after post-thinning statistics. */
public final class LesionPostThinningStandard13MembershipReferenceHarness {
    static final int N=17,C=13,NC=2,LESION=10;
    static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    static final float[][] MU=new float[NC][C],SD=new float[NC][C];
    static final float[] LC=new float[NC],LD=new float[NC],W=new float[NC],SP=new float[C],SM=new float[C];
    static final float[][][][] IMAGE=new float[N][N][N][NC],SHAPE=new float[N][N][N][C],MEM=new float[N][N][N][C],ENERGY=new float[N][N][N][C];
    static final boolean[][][][] DIRECT=new boolean[N][N][N][C];
    static final float[][][] LM=new float[N][N][N],BH=new float[N][N][N],OUT=new float[N][N][N],LREL=new float[N][N][N],DEN=new float[N][N][N],LE=new float[N][N][N],BHE=new float[N][N][N];
    static final int[][][] VENT=new int[N][N][N],GM=new int[N][N][N],BRAIN=new int[N][N][N];
    static final float[][] CLASS_SIM=new float[C][C]; static final float[] LESION_SIM=new float[C];
    static final float PRIOR=2.0f/13.0f,PRIOR_SCALE=0.1f,SMOOTH=(0.2f/6.0f)/14.0f,OUTLIER=0.1f,ZERO=1e-10f,INF=1e10f,EXP=1.2f;
    static final PowerTable POWER=new PowerTable(0f,4f,0.000001f,EXP);
    static void read(DataInputStream in)throws Exception{
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)MU[c][k]=in.readFloat();
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)SD[c][k]=in.readFloat();
        for(int c=0;c<NC;c++)LC[c]=in.readFloat();for(int c=0;c<NC;c++)LD[c]=in.readFloat();for(int c=0;c<NC;c++)W[c]=in.readFloat();
        for(int k=0;k<C;k++)SP[k]=in.readFloat();for(int k=0;k<C;k++)SM[k]=in.readFloat();
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            for(int c=0;c<NC;c++)IMAGE[x][y][z][c]=in.readFloat();
            for(int k=0;k<C;k++)SHAPE[x][y][z][k]=in.readFloat();
            for(int k=0;k<C;k++)DIRECT[x][y][z][k]=in.readInt()!=0;
            for(int k=0;k<C;k++)MEM[x][y][z][k]=in.readFloat();
            LM[x][y][z]=in.readFloat();BH[x][y][z]=in.readFloat();OUT[x][y][z]=in.readFloat();LREL[x][y][z]=in.readFloat();
            VENT[x][y][z]=in.readInt();GM[x][y][z]=in.readInt();BRAIN[x][y][z]=in.readInt();
        }
        if(in.read()!=-1)throw new IllegalArgumentException("trailing input bytes");
    }
    static int bits(float v){return Float.floatToRawIntBits(v);}
    static float sim(float d){float ratio=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;return 1.0f/(1.0f+ratio*d);}
    static float shapeFactor(float s){float twice=2.0f*s;return 1.0f<twice?1.0f:twice;}
    static void similarities(){
        for(int k=0;k<C;k++)for(int m=0;m<C;m++){
            if(k==m){CLASS_SIM[k][m]=0;continue;}float d=0;
            if(k==1||m==1){float q=MU[0][k]-MU[0][m],den=SD[0][k]+SD[0][m];d=(q*q/den)/den;}
            else for(int c=0;c<NC;c++){float q=MU[c][k]-MU[c][m],weighted=(W[c]*q)*q,den=SD[c][k]+SD[c][m];d+=(weighted/den)/den;}
            CLASS_SIM[k][m]=sim(d);
        }
        for(int k=0;k<C;k++){float d=0;for(int c=0;c<NC;c++){float q=LC[c]-MU[c][k],weighted=(W[c]*q)*q,den=LD[c]+SD[c][k];d+=(weighted/den)/den;}LESION_SIM[k]=k==LESION?0:sim(d);}
    }
    static void voxel(int x,int y,int z){
        float den=0;float[] faceL={LM[x+1][y][z],LM[x-1][y][z],LM[x][y+1][z],LM[x][y-1][z],LM[x][y][z+1],LM[x][y][z-1]};
        float[] faceB={BH[x+1][y][z],BH[x-1][y][z],BH[x][y+1][z],BH[x][y-1][z],BH[x][y][z+1],BH[x][y][z-1]};
        float[][] faceM={MEM[x+1][y][z],MEM[x-1][y][z],MEM[x][y+1][z],MEM[x][y-1][z],MEM[x][y][z+1],MEM[x][y][z-1]};
        for(int k=0;k<C;k++){
            float sh=SHAPE[x][y][z][k];if(!(sh>0)){MEM[x][y][z][k]=0;continue;}
            float raw,rawL=0;
            if(k==0){raw=0;for(int c=0;c<NC;c++)raw+=W[c]*(IMAGE[x][y][z][c]==0?0.0f:1.0f);}
            else if(k==1||k==2)raw=OUTLIER;
            else {raw=0;for(int c=0;c<NC;c++){float d=IMAGE[x][y][z][c]-MU[c][k],weighted=W[c]*d;raw+=weighted*d;}if(k==LESION){float d=IMAGE[x][y][z][1]-LC[1];raw=(float)raw;rawL=(float)POWER.lookup((double)(d*d),(double)EXP);}}
            float other=0,ls=0,bs=0,self=0;
            for(int f=0;f<6;f++){for(int m=0;m<C;m++)if(m!=k)other+=faceM[f][m]*faceM[f][m];ls+=faceL[f]*faceL[f];bs+=faceB[f]*faceB[f];}
            if(k==LESION)for(int f:new int[]{1,0,3,2,5,4})self+=faceM[f][k]*faceM[f][k];
            float sum=other+ls+bs,term=raw+(SM[k]*SMOOTH)*sum,energy=raw+((0.5f*SM[k])*SMOOTH)*sum,lt=rawL,le=0;
            if(k==LESION){float lsFor=other+self+bs;lt+=(SM[k]*SMOOTH)*lsFor;le=rawL+((0.5f*SM[k])*SMOOTH)*lsFor;}
            float ip=0,il=0;
            for(int m=0;m<C;m++){if(m!=k){float p=SP[m]*CLASS_SIM[k][m];p=p*SHAPE[x][y][z][m];p=p*SHAPE[x][y][z][m];ip+=p;}if(k==LESION){float p=SP[m]*LESION_SIM[m];p=p*SHAPE[x][y][z][m];p=p*SHAPE[x][y][z][m];il+=p;}}
            float p=SP[LESION]*LESION_SIM[k];p=p*SHAPE[x][y][z][LESION];p=p*SHAPE[x][y][z][LESION];ip+=p;
            term+=PRIOR*ip;energy+=((0.5f*PRIOR)*ip);if(k==LESION){lt+=PRIOR*il;le+=((0.5f*PRIOR)*il);}
            float cls=term>ZERO?1.0f/term:INF;if(!DIRECT[x][y][z][k])cls=shapeFactor(sh)*cls;MEM[x][y][z][k]=cls;ENERGY[x][y][z][k]=energy;
            if(k==LESION){float l=lt>ZERO?1.0f/lt:INF;if(!DIRECT[x][y][z][k])l=shapeFactor(sh)*l;LM[x][y][z]=l;LREL[x][y][z]=DIRECT[x][y][z][k]?1.0f:shapeFactor(sh);LE[x][y][z]=le;
                int vd=VENT[x][y][z];if(vd!=-1){float d=vd==0?0.5f:vd;LM[x][y][z]=d*LM[x][y][z]/2.0f;LREL[x][y][z]*=vd==0?d/2.0f:vd/2;}
                else if(GM[x][y][z]!=-1){int gd=GM[x][y][z];float d=gd==0?0.5f:gd;LM[x][y][z]=d*LM[x][y][z]/3.0f;LREL[x][y][z]*=gd==0?d/3.0f:gd/3;}
                int bd=BRAIN[x][y][z];if(bd!=-1){LM[x][y][z]=bd*LM[x][y][z]/25.0f;LREL[x][y][z]*=bd/25;}}
            if(k==4&&VENT[x][y][z]!=-1)MEM[x][y][z][k]=MEM[x][y][z][k]*VENT[x][y][z]/2.0f;
            if(k==LESION)den+=LM[x][y][z];den+=MEM[x][y][z][k];
        }
        DEN[x][y][z]=den;if(den>0){for(int k=0;k<C;k++){MEM[x][y][z][k]/=den;if(k==LESION)LM[x][y][z]/=den;}}
    }
    static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Integer.toUnsignedLong(bits(a[i])));}return s.append(']').toString();}
    static String flat4(float[][][][] a){float[] v=new float[N*N*N*a[0][0][0].length];int n=0;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<a[0][0][0].length;k++)v[n++]=a[x][y][z][k];return vec(v);}
    static String flat3(float[][][] a){float[] v=new float[N*N*N];int n=0;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)v[n++]=a[x][y][z];return vec(v);}
    static void run() {
        similarities();float objective=0;float[][][] running=new float[N][N][N];
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            voxel(x,y,z);if(DEN[x][y][z]>0)for(int k=0;k<C;k++){float q=ENERGY[x][y][z][k]*MEM[x][y][z][k];objective+=q*MEM[x][y][z][k];}
            float q=LE[x][y][z]*LM[x][y][z];q=q*LM[x][y][z];objective+=q;float b=BHE[x][y][z]*BH[x][y][z];b=b*BH[x][y][z];objective+=b;running[x][y][z]=objective;
        }
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"gm_flag\":true,\"class_membership_bits\":"+flat4(MEM)+",\"class_energy_bits\":"+flat4(ENERGY)+",\"lesion_membership_bits\":"+flat3(LM)+",\"black_hole_membership_bits\":"+flat3(BH)+",\"outlier_membership_bits\":"+flat3(OUT)+",\"lesion_relation_bits\":"+flat3(LREL)+",\"denominator_bits\":"+flat3(DEN)+",\"lesion_energy_bits\":"+flat3(LE)+",\"black_hole_energy_bits\":"+flat3(BHE)+",\"running_objective_bits\":"+flat3(running)+",\"objective_bits\":"+Integer.toUnsignedLong(bits(objective))+"}");
    }
    public static void main(String[] args)throws Exception{if(args.length!=1)throw new IllegalArgumentException("input binary required");try(DataInputStream in=new DataInputStream(new FileInputStream(args[0]))){read(in);}run();}
}
