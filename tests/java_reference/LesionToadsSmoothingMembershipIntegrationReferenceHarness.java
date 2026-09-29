/** Selected source formula composition: local smoothing through one voxel update. */
public final class LesionToadsSmoothingMembershipIntegrationReferenceHarness {
    private static final float[][] N={{0.10f,0.25f},{0.15f,0.30f},{0.20f,0.35f},{0.22f,0.37f},{0.12f,0.27f},{0.18f,0.33f}};
    private static final float[] L={0.10f,0.20f,0.30f,0.40f,0.50f,0.60f};
    private static final float[] B={0.05f,0.15f,0.25f,0.35f,0.45f,0.55f};
    private static final float ZERO=1.0e-10f,INF=1.0e10f;
    private static final class C {
        final String target,rel;final float[] term,shape,classEnergy;final float lesionTerm,bhTerm,lesionEnergy,bhEnergy;
        final float prior,smoothing,spread;final int vent,gm,bstem,maxVent,maxGM,maxBstem;final boolean gmFlag;
        C(String t,String r,float[] d,float[] s,float[] e,float ld,float bd,float le,float be,float p,float sm,
          int v,int g,int b,int mv,int mg,int mb,boolean gf,float sp){
            target=t;rel=r;term=d;shape=s;classEnergy=e;lesionTerm=ld;bhTerm=bd;lesionEnergy=le;bhEnergy=be;
            prior=p;smoothing=sm;vent=v;gm=g;bstem=b;maxVent=mv;maxGM=mg;maxBstem=mb;gmFlag=gf;spread=sp;
        }
    }
    private static final C[] CASES={
        new C("LESION","direct",new float[]{2.0f,3.0f},new float[]{0.45f,0.30f},new float[]{1.25f,1.75f},4.0f,0.0f,2.25f,0.0f,0.4f,0.6f,2,0,3,4,5,5,true,0.5f),
        new C("BLACKHOLE","prior",new float[]{1.25f,2.75f},new float[]{0.40f,0.20f},new float[]{0.75f,1.25f},3.5f,5.0f,0.5f,0.75f,0.8f,0.25f,-1,0,1,9,5,4,true,0.375f)
    };
    private static long bits(float v){return Float.floatToIntBits(v)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String matrix(float[][] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(vec(a[i]));}return s.append(']').toString();}
    private static float sf(float x){float twice=2.0f*x;return 1.0f<twice?1.0f:twice;}
    private static float inv(float e,String relation,float shape){float v=e>ZERO?1.0f/e:INF;return relation.equals("prior")?sf(shape)*v:v;}
    private static float[] smooth(C c,int k,String type,float lesionTerm,float bhTerm,float classEnergy,
        float lesionEnergy,float bhEnergy){
        float ngb=0,ngbL=0,ngbB=0,ngbS=0;
        for(int f=0;f<6;f++){
            for(int m=0;m<2;m++)if(m!=k)ngb+=N[f][m]*N[f][m];
            ngbL+=L[f]*L[f];ngbB+=B[f]*B[f];
        }
        if(type.equals("LESION")||type.equals("BLACKHOLE"))for(int f:new int[]{1,0,3,2,5,4})ngbS+=N[f][k]*N[f][k];
        float common=ngb+ngbL+ngbB;
        float data=c.term[k]+c.prior*c.smoothing*common;
        float energy=classEnergy+0.5f*c.prior*c.smoothing*common;
        float lt=lesionTerm,le=lesionEnergy,bt=bhTerm,be=bhEnergy;
        if(type.equals("LESION")||type.equals("BLACKHOLE")){
            float sum=ngb+ngbS+ngbB;
            lt+=c.prior*c.smoothing*sum;le+=0.5f*c.prior*c.smoothing*sum;
        }
        if(type.equals("BLACKHOLE")){
            float sum=ngb+ngbS+ngbL;
            bt+=c.prior*c.smoothing*sum;be+=0.5f*c.prior*c.smoothing*sum;
        }
        return new float[]{data,energy,lt,le,bt,be};
    }
    private static void run(C c,boolean first){
        String[] types={"OBJECT",c.target};float[] mem=new float[2],eng=new float[2];
        float lesion=0,bh=0,relation=0,lesionEnergy=0,bhEnergy=0;
        float[] ordinary=smooth(c,0,"OBJECT",0,0,c.classEnergy[0],0,0);
        mem[0]=inv(ordinary[0],"direct",c.shape[0]);eng[0]=ordinary[1];
        float[] target=smooth(c,1,c.target,c.lesionTerm,c.bhTerm,c.classEnergy[1],c.lesionEnergy,c.bhEnergy);
        mem[1]=inv(target[0],c.rel,c.shape[1]);eng[1]=target[1];
        lesion=inv(target[2],c.rel,c.shape[1]);lesionEnergy=target[3];
        if(c.target.equals("BLACKHOLE")){bh=inv(target[4],c.rel,c.shape[1]);bhEnergy=target[5];}
        relation=c.rel.equals("direct")?1.0f:sf(c.shape[1]);
        if(c.target.equals("LESION")){
            if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;relation*=n/c.maxVent;}
            else if(c.gmFlag&&c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;relation*=n/c.maxGM;}
            if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;relation*=(float)c.bstem/c.maxBstem;}
        }else{
            if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;relation*=n/c.maxVent;}
            else if(c.gmFlag&&c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;bh=n*bh/c.maxGM;relation*=n/c.maxGM;}
            if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;bh=c.bstem*bh/c.maxBstem;relation*=(float)c.bstem/c.maxBstem;}
        }
        float den=0.0f;for(int k=0;k<2;k++){if(types[k].equals("LESION"))den+=lesion;else if(types[k].equals("BLACKHOLE"))den+=lesion+bh;den+=mem[k];}
        float energy=0.0f;if(den>0.0f){for(int k=0;k<2;k++){mem[k]=mem[k]/den;if(types[k].equals("LESION")){lesion=lesion/den;}else if(types[k].equals("BLACKHOLE")){lesion=lesion/den;bh=bh/den;}energy+=eng[k]*mem[k]*mem[k];}}
        energy+=lesionEnergy*lesion*lesion+bhEnergy*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"target_type\":\""+c.target+"\",\"relation_type\":\""+c.rel+"\",\"class_term_bits\":"+vec(c.term)
            +",\"shape_bits\":"+vec(c.shape)+",\"class_energy_input_bits\":"+vec(c.classEnergy)
            +",\"lesion_term_bits\":"+bits(c.lesionTerm)+",\"black_hole_term_bits\":"+bits(c.bhTerm)
            +",\"lesion_energy_input_bits\":"+bits(c.lesionEnergy)+",\"black_hole_energy_input_bits\":"+bits(c.bhEnergy)
            +",\"smoothing_prior_bits\":"+bits(c.prior)+",\"smoothing_bits\":"+bits(c.smoothing)
            +",\"ventricle_distance\":"+c.vent+",\"gray_matter_distance\":"+c.gm+",\"brainstem_distance\":"+c.bstem
            +",\"max_ventricle_distance\":"+c.maxVent+",\"max_gray_matter_distance\":"+c.maxGM+",\"max_brainstem_distance\":"+c.maxBstem
            +",\"gm_flag\":"+c.gmFlag+",\"spread_bits\":"+bits(c.spread)+",\"neighbor_class_bits\":"+matrix(N)
            +",\"neighbor_lesion_bits\":"+vec(L)+",\"neighbor_black_hole_bits\":"+vec(B)
            +",\"class_membership_bits\":"+vec(mem)+",\"class_energy_bits\":"+vec(eng)
            +",\"lesion_membership_bits\":"+bits(lesion)+",\"black_hole_membership_bits\":"+bits(bh)
            +",\"lesion_relation_bits\":"+bits(relation)+",\"denominator_bits\":"+bits(den)+",\"energy_bits\":"+bits(energy)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships smoothing-to-membership composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);System.out.println("]}");
    }
}
