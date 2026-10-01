/** Selected source formula composition from six-face smoothing through normalization. */
public final class LesionToadsSmoothingPriorMembershipReferenceHarness {
    private static final float[][] NEIGHBORS={{0.10f,0.25f},{0.15f,0.30f},{0.20f,0.35f},{0.22f,0.37f},{0.12f,0.27f},{0.18f,0.33f}};
    private static final float[] NLESION={0.10f,0.20f,0.30f,0.40f,0.50f,0.60f};
    private static final float[] NBH={0.05f,0.15f,0.25f,0.35f,0.45f,0.55f};
    private static final float ZERO=1.0e-10f,INF=1.0e10f;
    private static final class C {
        final String target,relation;final boolean hasLesion;
        final float[] shapePriors,lesionSimilarity,blackHoleSimilarity,shapes,classTerms,classEnergies;
        final float[][] similarities;final float lesionTerm,bhTerm,lesionEnergy,bhEnergy,prior,smoothing,priorCoefficient,spread;
        final int vent,gm,bstem,maxVent,maxGM,maxBstem;
        C(String t,String r,boolean hl,float[] sp,float[][] sim,float[] ls,float[] bs,float[] sh,float[] ct,float[] ce,
          float lt,float bt,float le,float be,float p,float sm,float pc,int v,int g,int b,int mv,int mg,int mb,float spr){
            target=t;relation=r;hasLesion=hl;shapePriors=sp;similarities=sim;lesionSimilarity=ls;blackHoleSimilarity=bs;
            shapes=sh;classTerms=ct;classEnergies=ce;lesionTerm=lt;bhTerm=bt;lesionEnergy=le;bhEnergy=be;
            prior=p;smoothing=sm;priorCoefficient=pc;vent=v;gm=g;bstem=b;maxVent=mv;maxGM=mg;maxBstem=mb;spread=spr;
        }
    }
    private static final C[] CASES={
        new C("LESION","direct",true,new float[]{0.4f,0.6f},new float[][]{{0.9f,0.4f},{0.3f,0.8f}},new float[]{0.2f,0.7f},new float[]{0.6f,0.1f},new float[]{0.45f,0.3f},new float[]{1.5f,2.5f},new float[]{0.2f,0.3f},4.0f,0.0f,0.5f,0.0f,0.4f,0.6f,0.35f,2,0,3,4,5,5,0.5f),
        new C("BLACKHOLE","prior",true,new float[]{0.25f,0.75f},new float[][]{{0.7f,0.2f},{0.4f,0.9f}},new float[]{0.3f,0.8f},new float[]{0.9f,0.2f},new float[]{0.4f,0.2f},new float[]{1.25f,2.75f},new float[]{0.15f,0.25f},3.5f,5.0f,0.3f,0.4f,0.8f,0.25f,0.45f,-1,0,1,9,5,4,0.375f)
    };
    private static long bits(float x){return Float.floatToIntBits(x)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String matrix(float[][] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(vec(a[i]));}return s.append(']').toString();}
    private static float sf(float x){float twice=2.0f*x;return 1.0f<twice?1.0f:twice;}
    private static float inverse(float e,String relation,float shape){float v=e>ZERO?1.0f/e:INF;return relation.equals("prior")?sf(shape)*v:v;}
    private static float[] update(C c,int k){
        String type=k==0?"OBJECT":c.target,rel=k==0?"direct":c.relation;
        float ngb=0,ngbL=0,ngbB=0,ngbS=0;
        for(int f=0;f<6;f++){
            for(int m=0;m<2;m++)if(m!=k)ngb+=NEIGHBORS[f][m]*NEIGHBORS[f][m];
            ngbL+=NLESION[f]*NLESION[f];ngbB+=NBH[f]*NBH[f];
        }
        if(type.equals("LESION")||type.equals("BLACKHOLE"))for(int f:new int[]{1,0,3,2,5,4})ngbS+=NEIGHBORS[f][k]*NEIGHBORS[f][k];
        float common=ngb+ngbL+ngbB;
        float classTerm=c.classTerms[k]+c.prior*c.smoothing*common;
        float classEnergy=c.classEnergies[k]+0.5f*c.prior*c.smoothing*common;
        float lesionTerm=k==1?c.lesionTerm:0.0f,bhTerm=k==1?c.bhTerm:0.0f;
        float lesionEnergy=k==1?c.lesionEnergy:0.0f,bhEnergy=k==1?c.bhEnergy:0.0f;
        if(type.equals("LESION")||type.equals("BLACKHOLE")){
            float lsum=ngb+ngbS+ngbB;
            lesionTerm+=c.prior*c.smoothing*lsum;lesionEnergy+=0.5f*c.prior*c.smoothing*lsum;
        }
        if(type.equals("BLACKHOLE")){
            float bsum=ngb+ngbS+ngbL;
            bhTerm+=c.prior*c.smoothing*bsum;bhEnergy+=0.5f*c.prior*c.smoothing*bsum;
        }
        float ip=0,il=0,ib=0;
        for(int m=0;m<2;m++){
            float shape=c.shapes[m];
            if(m!=k)ip+=c.shapePriors[m]*c.similarities[k][m]*shape*shape;
            if(type.equals("LESION")||type.equals("BLACKHOLE"))il+=c.shapePriors[m]*c.lesionSimilarity[m]*shape*shape;
            if(type.equals("BLACKHOLE"))ib+=c.shapePriors[m]*c.blackHoleSimilarity[m]*shape*shape;
        }
        if(c.hasLesion){
            float shape=c.shapes[1];ip+=c.shapePriors[1]*c.lesionSimilarity[k]*shape*shape;
            if(c.target.equals("BLACKHOLE"))ip+=c.shapePriors[1]*c.blackHoleSimilarity[k]*shape*shape;
        }
        classTerm+=c.priorCoefficient*ip;classEnergy+=0.5f*c.priorCoefficient*ip;
        if(type.equals("LESION")||type.equals("BLACKHOLE")){
            lesionTerm+=c.priorCoefficient*il;lesionEnergy+=0.5f*c.priorCoefficient*il;
        }
        if(type.equals("BLACKHOLE")){
            bhTerm+=c.priorCoefficient*ib;bhEnergy+=0.5f*c.priorCoefficient*ib;
        }
        return new float[]{inverse(classTerm,rel,c.shapes[k]),classEnergy,lesionTerm,lesionEnergy,bhTerm,bhEnergy};
    }
    private static void run(C c,boolean first){
        String[] types={"OBJECT",c.target};float[] mem=new float[2],energy=new float[2];float lesion=0,bh=0,relation=0,le=0,be=0;
        for(int k=0;k<2;k++){
            float[] v=update(c,k);mem[k]=v[0];energy[k]=v[1];
            if(k==1){
                lesion=inverse(v[2],c.relation,c.shapes[k]);le=v[3];
                if(c.target.equals("BLACKHOLE")){bh=inverse(v[4],c.relation,c.shapes[k]);be=v[5];}
                relation=c.relation.equals("direct")?1.0f:sf(c.shapes[k]);
                if(c.target.equals("LESION")){
                    if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;if(c.vent==0)relation*=n/c.maxVent;else relation*=c.vent/c.maxVent;}
                    else if(c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;if(c.gm==0)relation*=n/c.maxGM;else relation*=c.gm/c.maxGM;}
                    if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;relation*=c.bstem/c.maxBstem;}
                }else{
                    if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;if(c.vent==0)relation*=n/c.maxVent;else relation*=c.vent/c.maxVent;}
                    else if(c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;bh=n*bh/c.maxGM;if(c.gm==0)relation*=n/c.maxGM;else relation*=c.gm/c.maxGM;}
                    if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;bh=c.bstem*bh/c.maxBstem;relation*=c.bstem/c.maxBstem;}
                }
            }
        }
        float den=0;for(int k=0;k<2;k++){if(types[k].equals("LESION"))den+=lesion;else if(types[k].equals("BLACKHOLE"))den+=lesion+bh;den+=mem[k];}
        float dist=0;if(den>0)for(int k=0;k<2;k++){mem[k]/=den;if(types[k].equals("LESION"))lesion/=den;else if(types[k].equals("BLACKHOLE")){lesion/=den;bh/=den;}dist+=energy[k]*mem[k]*mem[k];}
        dist+=le*lesion*lesion+be*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"target_type\":\""+c.target+"\",\"relation_type\":\""+c.relation+"\",\"has_lesion\":"+c.hasLesion
            +",\"shape_priors_bits\":"+vec(c.shapePriors)+",\"similarities_bits\":"+matrix(c.similarities)
            +",\"lesion_similarities_bits\":"+vec(c.lesionSimilarity)+",\"black_hole_similarities_bits\":"+vec(c.blackHoleSimilarity)
            +",\"shape_values_bits\":"+vec(c.shapes)+",\"class_terms_bits\":"+vec(c.classTerms)+",\"class_energies_bits\":"+vec(c.classEnergies)
            +",\"lesion_term_bits\":"+bits(c.lesionTerm)+",\"black_hole_term_bits\":"+bits(c.bhTerm)
            +",\"lesion_energy_input_bits\":"+bits(c.lesionEnergy)+",\"black_hole_energy_input_bits\":"+bits(c.bhEnergy)
            +",\"prior_bits\":"+bits(c.prior)+",\"smoothing_bits\":"+bits(c.smoothing)+",\"prior_coefficient_bits\":"+bits(c.priorCoefficient)
            +",\"ventricle_distance\":"+c.vent+",\"gray_matter_distance\":"+c.gm+",\"brainstem_distance\":"+c.bstem
            +",\"max_ventricle_distance\":"+c.maxVent+",\"max_gray_matter_distance\":"+c.maxGM+",\"max_brainstem_distance\":"+c.maxBstem
            +",\"spread_bits\":"+bits(c.spread)+",\"neighbor_class_bits\":"+matrix(NEIGHBORS)+",\"neighbor_lesion_bits\":"+vec(NLESION)
            +",\"neighbor_black_hole_bits\":"+vec(NBH)+",\"class_memberships_bits\":"+vec(mem)+",\"class_energy_bits\":"+vec(energy)
            +",\"lesion_membership_bits\":"+bits(lesion)+",\"black_hole_membership_bits\":"+bits(bh)+",\"lesion_relation_bits\":"+bits(relation)
            +",\"denominator_bits\":"+bits(den)+",\"energy_bits\":"+bits(dist)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships smoothing/prior composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);System.out.println("]}");
    }
}
