/** Selected one-voxel priorCoefficient source excerpt. */
public final class LesionToadsPriorEnergyReferenceHarness {
    private static final class C {
        final String type,lesionType;final int k,lesionClass;final boolean hasLesion;
        final float[] priors,similarity,lesionSimilarity,blackHoleSimilarity,shapes;
        final float coefficient,classTerm,lesionTerm,bhTerm,outmaskTerm,classEnergy,lesionEnergy,bhEnergy;
        C(String t,int i,String lt,int li,boolean hl,float[] p,float[] s,float[] ls,float[] bs,float[] sh,
          float pc,float ct,float l,float b,float o,float ce,float le,float be){
            type=t;k=i;lesionType=lt;lesionClass=li;hasLesion=hl;priors=p;similarity=s;
            lesionSimilarity=ls;blackHoleSimilarity=bs;shapes=sh;coefficient=pc;
            classTerm=ct;lesionTerm=l;bhTerm=b;outmaskTerm=o;classEnergy=ce;lesionEnergy=le;bhEnergy=be;
        }
    }
    private static final C[] CASES={
        new C("OBJECT",0,"LESION",2,true,new float[]{0.2f,0.3f,0.5f},new float[]{0.9f,0.4f,0.8f},new float[]{0.2f,0.5f,0.7f},new float[]{0.6f,0.3f,0.1f},new float[]{0.1f,0.25f,0.4f},0.35f,0.75f,0.25f,0.125f,0.0625f,0.5f,0.3f,0.2f),
        new C("LESION",1,"BLACKHOLE",2,true,new float[]{0.125f,0.25f,0.625f},new float[]{0.2f,0.7f,0.4f},new float[]{0.3f,0.8f,0.6f},new float[]{0.9f,0.2f,0.5f},new float[]{0.15f,0.35f,0.45f},0.6f,0.3f,0.4f,0.25f,0.125f,0.7f,0.8f,0.9f),
        new C("BLACKHOLE",2,"LESION",0,true,new float[]{0.4f,0.35f,0.25f},new float[]{0.3f,0.6f,0.2f},new float[]{0.7f,0.1f,0.9f},new float[]{0.2f,0.8f,0.4f},new float[]{0.12f,0.27f,0.38f},0.2f,0.15f,0.45f,0.55f,0.125f,0.1f,0.2f,0.3f),
        new C("OBJECT",1,"OBJECT",0,false,new float[]{0.0f,0.0f,0.0f},new float[]{0.2f,0.4f,0.6f},new float[]{0.3f,0.5f,0.7f},new float[]{0.1f,0.2f,0.3f},new float[]{0.0f,0.0f,0.0f},0.0f,0.25f,0.5f,0.75f,0.875f,1.0f,2.0f,3.0f),
        new C("OPTIMIZED_LESION",0,"BLACKHOLE",1,false,new float[]{0.1234567f,0.2345678f,0.6419755f},new float[]{0.8765432f,0.31415927f,0.27182818f},new float[]{0.14142136f,0.17320508f,0.2236068f},new float[]{0.24494898f,0.26457512f,0.28284273f},new float[]{0.11111111f,0.22222222f,0.33333334f},0.73123455f,0.91234565f,0.8123456f,0.7123456f,0.6123456f,0.5123457f,0.41234568f,0.31234568f),
        new C("OUTMASK",2,"OBJECT",0,false,new float[]{0.2f,0.3f,0.5f},new float[]{0.9f,0.4f,0.8f},new float[]{0.2f,0.5f,0.7f},new float[]{0.6f,0.3f,0.1f},new float[]{0.1f,0.25f,0.4f},0.35f,0.75f,0.25f,0.125f,0.5f,0.35f,0.3f,0.2f)
    };
    private static long bits(float v){return Float.floatToIntBits(v)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static void run(C c,boolean first){
        float ip=0,il=0,ib=0;
        for(int m=0;m<c.priors.length;m++){
            float shape=c.shapes[m];
            if(m!=c.k)ip+=c.priors[m]*c.similarity[m]*shape*shape;
            if(c.type.equals("LESION")||c.type.equals("OPTIMIZED_LESION")||c.type.equals("BLACKHOLE"))il+=c.priors[m]*c.lesionSimilarity[m]*shape*shape;
            if(c.type.equals("BLACKHOLE"))ib+=c.priors[m]*c.blackHoleSimilarity[m]*shape*shape;
        }
        if(c.hasLesion){
            float shape=c.shapes[c.lesionClass];
            ip+=c.priors[c.lesionClass]*c.lesionSimilarity[c.k]*shape*shape;
            if(c.lesionType.equals("BLACKHOLE"))ip+=c.priors[c.lesionClass]*c.blackHoleSimilarity[c.k]*shape*shape;
        }
        float classTerm=c.classTerm+c.coefficient*ip;
        float lesionTerm=c.lesionTerm,bhTerm=c.bhTerm;
        float outmaskTerm=c.outmaskTerm;
        if(c.type.equals("OUTMASK"))outmaskTerm+=c.coefficient*ip;
        if(c.type.equals("LESION")||c.type.equals("OPTIMIZED_LESION")||c.type.equals("BLACKHOLE"))lesionTerm+=c.coefficient*il;
        if(c.type.equals("BLACKHOLE"))bhTerm+=c.coefficient*ib;
        float ce=c.classEnergy+0.5f*c.coefficient*ip;
        float le=c.lesionEnergy,be=c.bhEnergy;
        if(c.type.equals("LESION")||c.type.equals("OPTIMIZED_LESION")||c.type.equals("BLACKHOLE"))le+=0.5f*c.coefficient*il;
        if(c.type.equals("BLACKHOLE"))be+=0.5f*c.coefficient*ib;
        if(!first)System.out.print(',');
        System.out.print("{\"class_type\":\""+c.type+"\",\"class_index\":"+c.k+",\"lesion_class_type\":\""+c.lesionType
            +"\",\"lesion_class_index\":"+c.lesionClass+",\"has_lesion\":"+c.hasLesion+",\"shape_priors_bits\":"+vec(c.priors)
            +",\"class_similarities_bits\":"+vec(c.similarity)+",\"lesion_similarities_bits\":"+vec(c.lesionSimilarity)
            +",\"black_hole_similarities_bits\":"+vec(c.blackHoleSimilarity)+",\"shape_values_bits\":"+vec(c.shapes)
            +",\"prior_coefficient_bits\":"+bits(c.coefficient)+",\"initial_class_term_bits\":"+bits(c.classTerm)
            +",\"initial_lesion_term_bits\":"+bits(c.lesionTerm)+",\"initial_black_hole_term_bits\":"+bits(c.bhTerm)
            +",\"initial_outmask_term_bits\":"+bits(c.outmaskTerm)
            +",\"initial_class_energy_bits\":"+bits(c.classEnergy)+",\"initial_lesion_energy_bits\":"+bits(c.lesionEnergy)
            +",\"initial_black_hole_energy_bits\":"+bits(c.bhEnergy)+",\"class_prior_bits\":"+bits(ip)
            +",\"lesion_prior_bits\":"+bits(il)+",\"black_hole_prior_bits\":"+bits(ib)
            +",\"class_term_bits\":"+bits(classTerm)+",\"lesion_term_bits\":"+bits(lesionTerm)+",\"black_hole_term_bits\":"+bits(bhTerm)
            +",\"outmask_term_bits\":"+bits(outmaskTerm)
            +",\"class_energy_bits\":"+bits(ce)+",\"lesion_energy_bits\":"+bits(le)+",\"black_hole_energy_bits\":"+bits(be)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships priorCoefficient excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);System.out.println("]}");
    }
}
