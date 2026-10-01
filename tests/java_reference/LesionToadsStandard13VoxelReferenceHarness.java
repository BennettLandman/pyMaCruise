/**
 * Fixed synthetic composition reference for VandyCRUISE
 * LesionToadSegmentation.computeMemberships(), source SHA-256
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * The pre-call state is the compact one-voxel harness contract in
 * tests/test_lesion_standard13_voxel.py; this harness mirrors the active
 * 13-class T1_MPRAGE+FLAIR branches and emits float bits for comparison.
 */
import edu.jhmi.rad.medic.utilities.PowerTable;
public final class LesionToadsStandard13VoxelReferenceHarness {
    static final int C=13;
    static final String[] N={"Background","Sulcal-CSF","Ventricles","Cerebellum-GM","Cerebrum-GM","Caudate","Thalamus","Putamen","Brainstem","Cerebellum-WM","OuterVentricular-WM","InterVentricular-WM","Frontal-WM"};
    static final String[] T={"MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","LESION","OBJECT","OBJECT"};
    static final float[] X={0.76f,0.9f}, LC={0.5f,0.65f}, LD={0.25f,0.3f}, W={1f,0.75f};
    static final float[][] MU=new float[2][C], SD=new float[2][C];
    static final float[] SH=new float[C], SP=new float[C], SM=new float[C], INIT=new float[C];
    static final boolean[] DIRECT={true,false,true,false,true,false,true,false,true,false,true,false,true};
    static final float[][] FACE=new float[6][C];
    static final float[] LF={0.05f,0.1f,0.15f,0.2f,0.25f,0.3f}, BF={0f,0.01f,0.02f,0.03f,0.04f,0.05f};
    static final float ZERO=1.0e-10f, INF=1.0e10f;
    static final float PRIOR=2.0f/13.0f, PRIOR_SCALE=0.1f, SMOOTH=(0.2f/6.0f)/14.0f;
    static final float OUTLIER=0.1f;
    static final float LESION_EXPONENT=1.2f;
    static final PowerTable LESION_POWER=new PowerTable(0.0f,4.0f,0.000001f,LESION_EXPONENT);

    static {
        for(int k=0;k<C;k++){
            MU[0][k]=0.2f+0.04f*k; MU[1][k]=0.3f+0.03f*k;
            SD[0][k]=0.3f+0.01f*k; SD[1][k]=0.4f+0.005f*k;
            SH[k]=0.4f+0.02f*k; SP[k]=0.5f+0.05f*k; SM[k]=0.8f+0.025f*k;
            INIT[k]=0.01f+0.001f*k;
            for(int f=0;f<6;f++) FACE[f][k]=0.02f+0.001f*(f*13+k);
        }
    }
    static int bits(float v){return Float.floatToIntBits(v);}
    static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Integer.toUnsignedLong(bits(a[i])));}return s.append(']').toString();}
    static String flatMatrix(float[][] a){float[] flat=new float[a.length*a[0].length];int n=0;for(float[] row:a)for(float v:row)flat[n++]=v;return vec(flat);}
    static float shapeFactor(float sh){float twice=2.0f*sh;return 1.0f<twice?1.0f:twice;}
    static float similarity(float distance){float ratio=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;return 1.0f/(1.0f+ratio*distance);}

    static void run(int ventDistance, int gmDistance, int brainstemDistance, boolean gmFlag){
        float[][] sim=new float[C][C]; float[] lsim=new float[C];
        for(int k=0;k<C;k++) for(int m=0;m<C;m++){
            if(k==m){sim[k][m]=0;continue;}
            float d=0;
            if(N[k].equals("Sulcal-CSF")||N[m].equals("Sulcal-CSF")){
                float delta=MU[0][k]-MU[0][m], den=SD[0][k]+SD[0][m]; d=(delta*delta/den)/den;
            }else for(int c=0;c<2;c++){
                float delta=MU[c][k]-MU[c][m], weighted=(W[c]*delta)*delta, den=SD[c][k]+SD[c][m];
                d+= (weighted/den)/den;
            }
            sim[k][m]=similarity(d);
        }
        for(int k=0;k<C;k++){
            float d=0;
            for(int c=0;c<2;c++){
                float delta=LC[c]-MU[c][k], weighted=(W[c]*delta)*delta, den=LD[c]+SD[c][k]; d+=(weighted/den)/den;
            }
            lsim[k]=k==10?0:similarity(d);
        }

        float[] mem=INIT.clone(), energy=new float[C];
        float lesion=0.125f,bh=0.0f,outlierMembership=0.0f,lesionRelation=0.0f;
        float lesionEnergy=0.0f,bhEnergy=0.0f,denominator=0.0f;
        for(int k=0;k<C;k++){
            float raw, rawL=0;
            if(k==0){raw=0;for(int c=0;c<2;c++)raw+=W[c]*(X[c]==0?0.0f:1.0f);}
            else if(k==1||k==2)raw=OUTLIER;
            else {
                raw=0;
                for(int c=0;c<2;c++){float delta=X[c]-MU[c][k], weighted=W[c]*delta;raw+=weighted*delta;}
                if(k==10){float delta=X[1]-LC[1];rawL=(float)LESION_POWER.lookup(delta*delta,LESION_EXPONENT);}
            }

            float other=0, ls=0, bs=0, self=0;
            for(int f=0;f<6;f++){
                for(int m=0;m<C;m++)if(m!=k)other+=FACE[f][m]*FACE[f][m];
                ls+=LF[f]*LF[f];bs+=BF[f]*BF[f];
            }
            if(k==10)for(int f:new int[]{1,0,3,2,5,4})self+=FACE[f][k]*FACE[f][k];
            float smSum=other+ls+bs;
            float clsTerm=raw+(SM[k]*SMOOTH)*smSum;
            float clsEnergy=raw+((0.5f*SM[k])*SMOOTH)*smSum;
            float lesTerm=rawL, lesEnergy=0;
            if(k==10){float aux=other+self+bs;lesTerm+=(SM[k]*SMOOTH)*aux;lesEnergy=rawL+((0.5f*SM[k])*SMOOTH)*aux;}

            float inv=0, invL=0;
            for(int m=0;m<C;m++){
                if(m!=k){float product=SP[m]*sim[k][m];product=product*SH[m];product=product*SH[m];inv+=product;}
                if(k==10){float product=SP[m]*lsim[m];product=product*SH[m];product=product*SH[m];invL+=product;}
            }
            float product=SP[10]*lsim[k];product=product*SH[10];product=product*SH[10];inv+=product;
            clsTerm+=PRIOR*inv; clsEnergy+=((0.5f*PRIOR)*inv);
            if(k==10){lesTerm+=PRIOR*invL;lesEnergy+=((0.5f*PRIOR)*invL);}


            float data=clsTerm>ZERO?1.0f/clsTerm:INF;
            if(!DIRECT[k])data=shapeFactor(SH[k])*data;
            mem[k]=data; energy[k]=clsEnergy;
            if(k==10){
                float les=lesTerm>ZERO?1.0f/lesTerm:INF;
                if(!DIRECT[k])les=shapeFactor(SH[k])*les;
                lesion=les;lesionRelation=DIRECT[k]?1.0f:shapeFactor(SH[k]);lesionEnergy=lesEnergy;
                // Match the selected source branch order: ventricle takes precedence;
                // otherwise GM is active only after the wrapper sets GMflag.
                if(ventDistance!=-1){float d=ventDistance==0?0.5f:ventDistance;lesion=d*lesion/2.0f;lesionRelation*=ventDistance==0?d/2.0f:ventDistance/2;}
                else if(gmFlag&&gmDistance!=-1){float d=gmDistance==0?0.5f:gmDistance;lesion=d*lesion/3.0f;lesionRelation*=gmDistance==0?d/3.0f:gmDistance/3;}
                if(brainstemDistance!=-1){lesion=brainstemDistance*lesion/25.0f;lesionRelation*=brainstemDistance/25;}
            }
            if(k==4&&ventDistance!=-1)mem[k]=mem[k]*1.0f/2.0f; // Cerebrum-GM ventricle-distance scaling
            if(k==10)denominator+=lesion;
            denominator+=mem[k];
        }
        if(denominator>0){for(int k=0;k<C;k++)mem[k]/=denominator;lesion/=denominator;}
        float objective=0;
        for(int k=0;k<C;k++){float q=energy[k]*mem[k];objective+=q*mem[k];}
        float q=lesionEnergy*lesion;objective+=q*lesion;
        q=bhEnergy*bh;objective+=q*bh;

        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"image_sha256\":\"8090f07543b5f45a9734f134a0cc927eb9d336f6e3c23b2b6b483a89322e2c40\"");
        System.out.print(",\"class_names\":[");for(int i=0;i<C;i++){if(i>0)System.out.print(',');System.out.print("\""+N[i]+"\"");}
        System.out.print("],\"class_types\":[");for(int i=0;i<C;i++){if(i>0)System.out.print(',');System.out.print("\""+T[i]+"\"");}
        System.out.print("],\"profile_bits\":{\"prior_coefficient\":"+Integer.toUnsignedLong(bits(PRIOR))+",\"prior_scale\":"+Integer.toUnsignedLong(bits(PRIOR_SCALE))+",\"class_smoothing\":"+Integer.toUnsignedLong(bits(SMOOTH))+"}");
        System.out.print(",\"distance_maxima\":{\"ventricle\":2,\"gray_matter\":3,\"brainstem\":25}");
        System.out.print(",\"class_membership_bits\":"+vec(mem)+",\"lesion_membership_bits\":"+Integer.toUnsignedLong(bits(lesion))
            +",\"black_hole_membership_bits\":"+Integer.toUnsignedLong(bits(bh)) +",\"outlier_membership_bits\":"+Integer.toUnsignedLong(bits(outlierMembership))
            +",\"lesion_relation_bits\":"+Integer.toUnsignedLong(bits(lesionRelation))+",\"class_energy_bits\":"+vec(energy)
            +",\"lesion_energy_bits\":"+Integer.toUnsignedLong(bits(lesionEnergy))+",\"black_hole_energy_bits\":"+Integer.toUnsignedLong(bits(bhEnergy))
            +",\"denominator_bits\":"+Integer.toUnsignedLong(bits(denominator))+",\"objective_bits\":"+Integer.toUnsignedLong(bits(objective))
            +",\"class_similarity_bits\":"+flatMatrix(sim)+",\"lesion_similarity_bits\":"+vec(lsim)+"}");
    }
    public static void main(String[] args){
        int vent=1, gm=-1, brainstem=10; boolean gmFlag=false;
        if(args.length!=0){if(args.length!=4)throw new IllegalArgumentException("usage: [ventDistance gmDistance brainstemDistance gmFlag]");vent=Integer.parseInt(args[0]);gm=Integer.parseInt(args[1]);brainstem=Integer.parseInt(args[2]);gmFlag=Boolean.parseBoolean(args[3]);}
        run(vent,gm,brainstem,gmFlag);
    }
}
