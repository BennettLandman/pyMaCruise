import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;
import edu.jhmi.rad.medic.utilities.PowerTable;

/**
 * Source-ordered selected VandyCRUISE standard-13 computeMemberships fixture.
 * Source SHA-256: 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * The fixed pre-call maps and RAS image formula are also recorded in
 * tests/test_lesion_standard13_volume.py. The 15 PowerTable outputs are
 * independently pinned from the selected FLAIR exponent-1.2 lookup.
 */
public final class LesionToadsInitialThinningCompositionReferenceHarness {
    static final int N=17,C=13,LESION_CLASS=10,EMPTY=-1;
    static final float VOLUME_RATIO=0.4f,LIMIT=1000.0f;
    static final String[] NAMES={"Background","Sulcal-CSF","Ventricles","Cerebellum-GM","Cerebrum-GM","Caudate","Thalamus","Putamen","Brainstem","Cerebellum-WM","OuterVentricular-WM","InterVentricular-WM","Frontal-WM"};
    static final float[] W={1.0f,0.75f}, LC={0.5f,0.65f}, LD={0.25f,0.3f};
    static final float[][] MU=new float[2][C], SD=new float[2][C];
    static final float[] SP=new float[C], SM=new float[C];
    static final float LESION_EXPONENT=1.2f;
    static final PowerTable LESION_POWER=new PowerTable(0.0f,4.0f,0.000001f,LESION_EXPONENT);
    static final float PRIOR=2.0f/13.0f, SCALE=0.1f, SMOOTH=(0.2f/6.0f)/14.0f, OUTLIER=0.1f, ZERO=1.0e-10f, INF=1.0e10f;
    static final float[][][][] IMAGE=new float[N][N][N][2], SHAPE=new float[N][N][N][C], MEM=new float[N][N][N][C], ENERGY=new float[N][N][N][C];
    static final boolean[][][][] DIRECT=new boolean[N][N][N][C];
    static final float[][][] LESION=new float[N][N][N], BH=new float[N][N][N], OUTMEM=new float[N][N][N], LREL=new float[N][N][N], LE=new float[N][N][N], BHE=new float[N][N][N], DEN=new float[N][N][N], RUNNING=new float[N][N][N];
    static final short[][][] SEG=new short[N][N][N], SKELETON=new short[N][N][N];
    static final boolean[][][][] AVAILABLE=new boolean[N][N][N][C];
    static final float[][][] ORDERING=new float[N][N][N];

    static {
        for(int k=0;k<C;k++){
            MU[0][k]=0.2f+0.04f*k; MU[1][k]=0.3f+0.03f*k;
            SD[0][k]=0.3f+0.01f*k; SD[1][k]=0.4f+0.005f*k;
            SP[k]=0.5f+0.05f*k; SM[k]=0.8f+0.025f*k;
        }
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int ti=(x*7+y*5+z*3+((x*y)%3)*2)%16, fi=(x*3+y*7+z*5+((y*z)%4)*3)%15;
            int fold=8+((3*x+5*y)%5);
            float folded=(z>=fold)?1.0f:0.0f;
            IMAGE[x][y][z][0]=0.4f+0.025f*ti+0.018f*folded;
            IMAGE[x][y][z][1]=0.55f+0.02f*fi+0.027f*folded;
            for(int k=0;k<C;k++){
                int offset=(x+2*y+3*z)%7;
                SHAPE[x][y][z][k]=(0.4f+0.02f*k)+0.001f*offset;
                DIRECT[x][y][z][k]=((x+2*y+3*z+k)%2)==0;
            }
        }
    }
    static int bits(float v){return Float.floatToIntBits(v);}
    static String vec(float[] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(float v:a){if(!first)s.append(',');first=false;s.append(Integer.toUnsignedLong(bits(v)));}return s.append(']').toString();}
    static String flat4(float[][][][] a){float[] f=new float[N*N*N*a[0][0][0].length];int i=0;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<a[0][0][0].length;k++)f[i++]=a[x][y][z][k];return vec(f);}
    static String flat3(float[][][] a){float[] f=new float[N*N*N];int i=0;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)f[i++]=a[x][y][z];return vec(f);}
    static float sf(float shape){float twice=2.0f*shape;return 1.0f<twice?1.0f:twice;}
    static float sim(float dist){float ratio=(1.0f-SCALE)/SCALE;return 1.0f/(1.0f+ratio*dist);}
    static float[][] CLASS_SIM=new float[C][C]; static float[] LESION_SIM=new float[C];

    static void similarities(){
        for(int k=0;k<C;k++)for(int m=0;m<C;m++){
            if(k==m){CLASS_SIM[k][m]=0;continue;}float dist=0;
            if(k==1||m==1){float d=MU[0][k]-MU[0][m],den=SD[0][k]+SD[0][m];dist=(d*d/den)/den;}
            else for(int c=0;c<2;c++){float d=MU[c][k]-MU[c][m],weighted=(W[c]*d)*d,den=SD[c][k]+SD[c][m];dist+=(weighted/den)/den;}
            CLASS_SIM[k][m]=sim(dist);
        }
        for(int k=0;k<C;k++){
            float dist=0;for(int c=0;c<2;c++){float d=LC[c]-MU[c][k],weighted=(W[c]*d)*d,den=LD[c]+SD[c][k];dist+=(weighted/den)/den;}
            LESION_SIM[k]=k==10?0:sim(dist);
        }
    }
    static int vent(int x,int y,int z){return ((x+y+z)%4==0)?-1:(x+2*y+z)%3;}
    static int gm(int x,int y,int z){return -1;}
    static int bstem(int x,int y,int z){return ((x+z)%5==0)?-1:(2*x+y+z)%26;}

    static void updateVoxel(int x,int y,int z){
        float[] faceL={LESION[x+1][y][z],LESION[x-1][y][z],LESION[x][y+1][z],LESION[x][y-1][z],LESION[x][y][z+1],LESION[x][y][z-1]};
        float[] faceB={BH[x+1][y][z],BH[x-1][y][z],BH[x][y+1][z],BH[x][y-1][z],BH[x][y][z+1],BH[x][y][z-1]};
        float[][] faceM={MEM[x+1][y][z],MEM[x-1][y][z],MEM[x][y+1][z],MEM[x][y-1][z],MEM[x][y][z+1],MEM[x][y][z-1]};
        float[] shape=SHAPE[x][y][z];float den=0;
        for(int k=0;k<C;k++){
            if(!(shape[k]>0)){MEM[x][y][z][k]=0;continue;}
            float raw,rawL=0;
            if(k==0){raw=0;for(int c=0;c<2;c++)raw+=W[c]*(IMAGE[x][y][z][c]==0?0.0f:1.0f);}
            else if(k==1||k==2)raw=OUTLIER;
            else {raw=0;for(int c=0;c<2;c++){float d=IMAGE[x][y][z][c]-MU[c][k],weighted=W[c]*d;raw+=weighted*d;}if(k==10){float d=IMAGE[x][y][z][1]-LC[1];rawL=(float)LESION_POWER.lookup(d*d,LESION_EXPONENT);}}
            float other=0,ls=0,bs=0,self=0;
            for(int f=0;f<6;f++){for(int m=0;m<C;m++)if(m!=k)other+=faceM[f][m]*faceM[f][m];ls+=faceL[f]*faceL[f];bs+=faceB[f]*faceB[f];}
            if(k==10)for(int f:new int[]{1,0,3,2,5,4})self+=faceM[f][k]*faceM[f][k];
            float sum=other+ls+bs;
            float term=raw+(SM[k]*SMOOTH)*sum;
            float energy=raw+((0.5f*SM[k])*SMOOTH)*sum;
            float lesionTerm=rawL, lesionEnergy=0;
            if(k==10){float lsFor=other+self+bs;lesionTerm+=(SM[k]*SMOOTH)*lsFor;lesionEnergy=rawL+((0.5f*SM[k])*SMOOTH)*lsFor;}
            float ip=0,il=0;
            for(int m=0;m<C;m++){
                if(m!=k){float p=SP[m]*CLASS_SIM[k][m];p=p*shape[m];p=p*shape[m];ip+=p;}
                if(k==10){float p=SP[m]*LESION_SIM[m];p=p*shape[m];p=p*shape[m];il+=p;}
            }
            float p=SP[10]*LESION_SIM[k];p=p*shape[10];p=p*shape[10];ip+=p;
            term+=PRIOR*ip;energy+=((0.5f*PRIOR)*ip);
            if(k==10){lesionTerm+=PRIOR*il;lesionEnergy+=((0.5f*PRIOR)*il);}
            float cls=term>ZERO?1.0f/term:INF;
            if(!DIRECT[x][y][z][k])cls=sf(shape[k])*cls;
            MEM[x][y][z][k]=cls;ENERGY[x][y][z][k]=energy;
            if(k==10){
                float l=lesionTerm>ZERO?1.0f/lesionTerm:INF;
                if(!DIRECT[x][y][z][k])l=sf(shape[k])*l;
                LESION[x][y][z]=l;LREL[x][y][z]=DIRECT[x][y][z][k]?1.0f:sf(shape[k]);LE[x][y][z]=lesionEnergy;
                int vd=vent(x,y,z);if(vd!=-1){float num=vd==0?0.5f:vd;LESION[x][y][z]=num*LESION[x][y][z]/2.0f;if(vd==0)LREL[x][y][z]*=num/2.0f;else LREL[x][y][z]*=vd/2;}
                int bd=bstem(x,y,z);if(bd!=-1){LESION[x][y][z]=bd*LESION[x][y][z]/25.0f;LREL[x][y][z]*=bd/25;}
            }
            if(k==4&&vent(x,y,z)!=-1)MEM[x][y][z][k]=MEM[x][y][z][k]*vent(x,y,z)/2.0f;
            if(k==10)den+=LESION[x][y][z];
            den+=MEM[x][y][z][k];
        }
        DEN[x][y][z]=den;
        if(den>0){
            for(int k=0;k<C;k++){MEM[x][y][z][k]/=den;if(k==10)LESION[x][y][z]/=den;}
        }
    }
    static void run(){
        similarities();float objective=0;
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            updateVoxel(x,y,z);
            if(DEN[x][y][z]>0)for(int k=0;k<C;k++){float q=ENERGY[x][y][z][k]*MEM[x][y][z][k];q=q*MEM[x][y][z][k];objective+=q;}
            float q=LE[x][y][z]*LESION[x][y][z];q=q*LESION[x][y][z];float b=BHE[x][y][z]*BH[x][y][z];b=b*BH[x][y][z];objective+=q+b;RUNNING[x][y][z]=objective;
        }
        runThinning();
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"image_sha256\":\"59d0a14b76ebe7f1fd045ef45add887b2fc8978015cdff5649d25cb313835f31\",\"wrapper_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"power_table_source_sha256\":\"5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678\",\"power_table_float_exponent_bits\":\"3f99999a\",\"heap_source_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"topology_lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"source_call_order\":[\"computeMemberships(GMflag=false)\",\"setGMflag\",\"Eprev=initial_energy\",\"Eratio=1.0f\",\"initThinningLabels\",\"propagateThinning\"]");
        System.out.print(",\"class_membership_bits\":"+flat4(MEM)+",\"class_energy_bits\":"+flat4(ENERGY)
            +",\"lesion_membership_bits\":"+flat3(LESION)+",\"black_hole_membership_bits\":"+flat3(BH)
            +",\"outlier_membership_bits\":"+flat3(OUTMEM)+",\"lesion_relation_bits\":"+flat3(LREL)
            +",\"denominator_bits\":"+flat3(DEN)+",\"lesion_energy_bits\":"+flat3(LE)
            +",\"black_hole_energy_bits\":"+flat3(BHE)+",\"running_objective_bits\":"+flat3(RUNNING)
            +",\"objective_bits\":"+Integer.toUnsignedLong(bits(objective))
            +",\"initial_volume_bits\":"+VOLUME_BITS+",\"minimum_volume_bits\":"+MINIMUM_BITS
            +",\"initial_queue\":"+INITIAL_QUEUE+",\"events\":"+EVENTS+",\"remaining_queue\":"+REMAINING
            +",\"final_segmentation\":"+flatShort(SEG)+",\"final_skeleton\":"+flatShort(SKELETON)
            +",\"final_available\":"+flatBool(AVAILABLE)+",\"final_ordering_bits\":"+flat3(ORDERING)
            +",\"final_volume_bits\":"+FINAL_VOLUME_BITS+"}");
    }
    static String VOLUME_BITS="[]",MINIMUM_BITS="[]",INITIAL_QUEUE="[]",EVENTS="[]",REMAINING="[]",FINAL_VOLUME_BITS="[]";
    static String flatShort(short[][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!first)s.append(',');first=false;s.append(a[x][y][z]);}return s.append(']').toString();}
    static String flatBool(boolean[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){if(!first)s.append(',');first=false;s.append(a[x][y][z][k]);}return s.append(']').toString();}
    static String encItem(float p,int x,int y,int z,int k){return "["+bits(p)+","+x+","+y+","+z+","+k+"]";}
    static int bestThinning(int x,int y,int z){int best=0;float running=0.0f;for(int k=0;k<C;k++){float score=k==LESION_CLASS?((MEM[x][y][z][k]+LESION[x][y][z])+BH[x][y][z]):MEM[x][y][z][k];if(score>running){best=k;running=MEM[x][y][z][k];}}return best;}
    static float ratioThinning(int x,int y,int z,int current,int best){if(current==best)return 1.0f;float cur=current==LESION_CLASS?((MEM[x][y][z][current]+LESION[x][y][z])+BH[x][y][z]):MEM[x][y][z][current];float max=best==LESION_CLASS?((MEM[x][y][z][best]+LESION[x][y][z])+BH[x][y][z]):MEM[x][y][z][best];return (max+ZERO)/(cur+ZERO);}
    static boolean boundary(int x,int y,int z,int k){return SEG[x-1][y][z]!=k||SEG[x+1][y][z]!=k||SEG[x][y-1][z]!=k||SEG[x][y+1][z]!=k||SEG[x][y][z-1]!=k||SEG[x][y][z+1]!=k;}
    static void runThinning(){
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){int fold=8+((3*x+5*y)%5);SEG[x][y][z]=(short)(z<fold?1:12);boolean border=x==0||x==N-1||y==0||y==N-1||z==0||z==N-1;SKELETON[x][y][z]=(short)(border?7:EMPTY);for(int k=0;k<C;k++)AVAILABLE[x][y][z][k]=((x+y+z+k)%2)==0;}
        float[] volume=new float[C];
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){int cur=SEG[x][y][z];SKELETON[x][y][z]=(short)cur;for(int k=0;k<C;k++){if(cur==k){if(k==LESION_CLASS)volume[k]+=(MEM[x][y][z][k]+LESION[x][y][z])+BH[x][y][z];else volume[k]+=MEM[x][y][z][k];}AVAILABLE[x][y][z][k]=cur==k;}}
        float[] initialVolume=volume.clone();BinaryHeap4D tree=new BinaryHeap4D(N*N*N*C,BinaryHeap4D.MAXTREE);List<String> inserted=new ArrayList<>();
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){int k=SEG[x][y][z];if(boundary(x,y,z,k)){AVAILABLE[x][y][z][k]=false;float priority=-1.0f/ratioThinning(x,y,z,k,bestThinning(x,y,z));tree.addValue(priority,x,y,z,k);inserted.add(encItem(priority,x,y,z,k));}}
        float[] minimum=new float[C];for(int k=0;k<C;k++)minimum[k]=Numerics.max(100.0f,VOLUME_RATIO*volume[k]);List<String> events=new ArrayList<>();float value=0.0f;
        while(tree.isNotEmpty()&&value>-LIMIT){short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();value=tree.getFirst();tree.removeFirst();boolean processed=false;int candidate=-1;String topo="null";List<String> added=new ArrayList<>();if(volume[k]>minimum[k]&&SKELETON[x][y][z]==k){int id=bestThinning(x,y,z);boolean direct=false;int next=-1;float bestPrior=0.0f;if(id!=k){for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;short neighbor=SEG[xi][yj][zl];if(neighbor==id)direct=true;else if(neighbor!=k&&SHAPE[xi][yj][zl][id]>bestPrior){next=neighbor;bestPrior=SHAPE[xi][yj][zl][id];}}if(direct)candidate=id;else if(next!=-1)candidate=next;if(candidate!=-1){boolean accepted=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(SEG,x,y,z,(short)candidate);topo=Boolean.toString(accepted);if(accepted)SEG[x][y][z]=(short)candidate;}}
            AVAILABLE[x][y][z][k]=false;ORDERING[x][y][z]=value;float dec=k==LESION_CLASS?((MEM[x][y][z][k]+LESION[x][y][z])+BH[x][y][z]):MEM[x][y][z][k];volume[k]-=dec;SKELETON[x][y][z]=(short)EMPTY;processed=true;for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(AVAILABLE[xi][yj][zl][k]){float priority=value-1.0f/ratioThinning(xi,yj,zl,k,bestThinning(xi,yj,zl));if(volume[k]>minimum[k]&&priority>-LIMIT){AVAILABLE[xi][yj][zl][k]=false;tree.addValue(priority,xi,yj,zl,k);added.add(encItem(priority,xi,yj,zl,k));}}}}
            events.add("["+bits(value)+","+x+","+y+","+z+","+k+","+processed+","+(candidate<0?"null":candidate)+","+topo+",["+String.join(",",added)+"]]");}
        List<String> remaining=new ArrayList<>();while(tree.isNotEmpty()){remaining.add(encItem(tree.getFirst(),tree.getFirstX(),tree.getFirstY(),tree.getFirstZ(),tree.getFirstK()));tree.removeFirst();}
        VOLUME_BITS=vec(initialVolume);MINIMUM_BITS=vec(minimum);FINAL_VOLUME_BITS=vec(volume);INITIAL_QUEUE="["+String.join(",",inserted)+"]";EVENTS="["+String.join(",",events)+"]";REMAINING="["+String.join(",",remaining)+"]";
    }
    public static void main(String[] args){run();}
}
