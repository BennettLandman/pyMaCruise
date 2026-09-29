/** Selected one-voxel six-face smoothing-energy source excerpt. */
public final class LesionToadsSmoothingEnergyReferenceHarness {
    private static final class C {
        final String type;final int index;final float[][] classes;final float[] lesion,bh;
        final float prior,smoothing,classTerm,lesionTerm,bhTerm,outmaskTerm,classEnergy,lesionEnergy,bhEnergy;
        C(String t,int k,float[][] m,float[] l,float[] b,float p,float s,float ct,float lt,float bt,float ot,float ce,float le,float be){
            type=t;index=k;classes=m;lesion=l;bh=b;prior=p;smoothing=s;classEnergy=ce;
            lesionEnergy=le;bhEnergy=be;classTerm=ct;lesionTerm=lt;bhTerm=bt;outmaskTerm=ot;
        }
    }
    private static final float[][] M={
        {0.125f,0.25f,0.375f},{0.2f,0.3f,0.4f},{0.45f,0.15f,0.05f},
        {0.11f,0.22f,0.33f},{0.17f,0.27f,0.37f},{0.19f,0.29f,0.39f}
    };
    private static final float[] L={0.13f,0.23f,0.31f,0.41f,0.51f,0.61f};
    private static final float[] B={0.07f,0.17f,0.27f,0.37f,0.47f,0.57f};
    private static final C[] CASES={
        new C("OBJECT",0,M,L,B,0.6f,0.7f,0.11f,0.21f,0.31f,0.41f,0.2f,0.3f,0.4f),
        new C("LESION",1,M,L,B,0.6f,0.7f,0.11f,0.21f,0.31f,0.41f,0.2f,0.3f,0.4f),
        new C("BLACKHOLE",2,M,L,B,0.75f,0.35f,0.125f,0.225f,0.325f,0.425f,0.25f,0.35f,0.45f),
        new C("OUTMASK",0,M,L,B,0.25f,0.8f,0.4f,0.3f,0.2f,0.1f,0.14f,0.24f,0.34f)
    };
    private static long bits(float v){return Float.floatToIntBits(v)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();}
    private static String matrix(float[][] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(vec(a[i]));}return b.append(']').toString();}
    private static void run(C c,boolean first){
        float ngb=0.0f,ngbL=0.0f,ngbB=0.0f,ngbS=0.0f;
        for(int f=0;f<6;f++){
            for(int m=0;m<c.classes[f].length;m++)if(m!=c.index)ngb+=c.classes[f][m]*c.classes[f][m];
            ngbL+=c.lesion[f]*c.lesion[f];ngbB+=c.bh[f]*c.bh[f];
        }
        if(c.type.equals("LESION")||c.type.equals("OPTIMIZED_LESION")||c.type.equals("BLACKHOLE")){
            int[] sourceOrder={1,0,3,2,5,4};
            for(int f:sourceOrder)ngbS+=c.classes[f][c.index]*c.classes[f][c.index];
        }
        float classSum=ngb+ngbL+ngbB;
        float classTerm=c.classTerm+c.prior*c.smoothing*classSum;
        float classEnergy=c.classEnergy+0.5f*c.prior*c.smoothing*classSum;
        float lesionTerm=c.lesionTerm,bhTerm=c.bhTerm,outmaskTerm=c.outmaskTerm;
        float lesionEnergy=c.lesionEnergy,bhEnergy=c.bhEnergy;
        if(c.type.equals("LESION")||c.type.equals("OPTIMIZED_LESION")||c.type.equals("BLACKHOLE")){
            float lsum=ngb+ngbS+ngbB;
            lesionTerm+=c.prior*c.smoothing*lsum;
            lesionEnergy+=0.5f*c.prior*c.smoothing*lsum;
        }
        if(c.type.equals("BLACKHOLE")){
            float bsum=ngb+ngbS+ngbL;
            bhTerm+=c.prior*c.smoothing*bsum;
            bhEnergy+=0.5f*c.prior*c.smoothing*bsum;
        }
        if(c.type.equals("OUTMASK"))outmaskTerm+=c.prior*c.smoothing*(ngb+ngbL);
        if(!first)System.out.print(',');
        System.out.print("{\"class_type\":\""+c.type+"\",\"class_index\":"+c.index+",\"class_neighbor_bits\":"+matrix(c.classes)
            +",\"lesion_neighbor_bits\":"+vec(c.lesion)+",\"black_hole_neighbor_bits\":"+vec(c.bh)
            +",\"smoothing_prior_bits\":"+bits(c.prior)+",\"smoothing_bits\":"+bits(c.smoothing)
            +",\"initial_class_data_term_bits\":"+bits(c.classTerm)+",\"initial_lesion_data_term_bits\":"+bits(c.lesionTerm)
            +",\"initial_black_hole_data_term_bits\":"+bits(c.bhTerm)+",\"initial_outmask_data_term_bits\":"+bits(c.outmaskTerm)
            +",\"initial_class_energy_bits\":"+bits(c.classEnergy)+",\"initial_lesion_energy_bits\":"+bits(c.lesionEnergy)
            +",\"initial_black_hole_energy_bits\":"+bits(c.bhEnergy)
            +",\"class_data_term_bits\":"+bits(classTerm)+",\"lesion_data_term_bits\":"+bits(lesionTerm)
            +",\"black_hole_data_term_bits\":"+bits(bhTerm)+",\"outmask_data_term_bits\":"+bits(outmaskTerm)
            +",\"class_energy_bits\":"+bits(classEnergy)+",\"lesion_energy_bits\":"+bits(lesionEnergy)
            +",\"black_hole_energy_bits\":"+bits(bhEnergy)
            +",\"other_class_neighbor_squares_bits\":"+bits(ngb)+",\"lesion_neighbor_squares_bits\":"+bits(ngbL)
            +",\"black_hole_neighbor_squares_bits\":"+bits(ngbB)+",\"current_class_neighbor_squares_bits\":"+bits(ngbS)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships six-face smoothing excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"face_order\":[\"X+\",\"X-\",\"Y+\",\"Y-\",\"Z+\",\"Z-\"],\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);
        System.out.println("]}");
    }
}
