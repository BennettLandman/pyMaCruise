/** Selected per-voxel LesionTOADS membership formula composition. */
public final class LesionToadsMembershipIntegrationReferenceHarness {
    private static final float ZERO=1.0e-10f,INF=1.0e10f;
    private static final class C {
        final String[] types,relations;final float[] energies,shapes,lesionE,blackHoleE;
        final int vent,gm,bstem,maxVent,maxGM,maxBstem;final boolean gmFlag;final float spread;
        C(String[] t,String[] r,float[] e,float[] s,float[] le,float[] be,
          int v,int g,int b,int mv,int mg,int mb,boolean gf,float sp){
            types=t;relations=r;energies=e;shapes=s;lesionE=le;blackHoleE=be;
            vent=v;gm=g;bstem=b;maxVent=mv;maxGM=mg;maxBstem=mb;gmFlag=gf;spread=sp;
        }
    }
    private static final C[] CASES={
        new C(new String[]{"OBJECT","LESION"},new String[]{"direct","direct"},new float[]{2.0f,3.0f},new float[]{0.5f,0.25f},new float[]{0.0f,4.0f},new float[]{0.0f,0.0f},2,3,4,4,5,5,true,0.5f),
        new C(new String[]{"OBJECT","LESION"},new String[]{"prior","prior"},new float[]{1.5f,2.0f},new float[]{0.25f,0.2f},new float[]{0.0f,4.0f},new float[]{0.0f,0.0f},-1,0,-1,9,5,5,true,0.375f),
        new C(new String[]{"OBJECT","BLACKHOLE"},new String[]{"direct","direct"},new float[]{2.5f,3.5f},new float[]{0.5f,0.3f},new float[]{0.0f,4.0f},new float[]{0.0f,5.0f},2,0,3,4,5,5,true,0.5f),
        new C(new String[]{"OBJECT","LESION","BLACKHOLE"},new String[]{"direct","prior","direct"},new float[]{1.25f,2.75f,3.125f},new float[]{0.5f,0.35f,0.4f},new float[]{0.0f,4.25f,2.5f},new float[]{0.0f,0.0f,6.0f},-1,2,1,4,5,4,true,0.25f),
        new C(new String[]{"LESION"},new String[]{"prior"},new float[]{1.0f},new float[]{0.0f},new float[]{2.0f},new float[]{0.0f},-1,-1,-1,0,0,0,false,0.5f)
    };
    private static long bits(float v){return Float.floatToIntBits(v)&0xffffffffL;}
    private static String fbits(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();}
    private static String strings(String[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('"').append(a[i]).append('"');}return b.append(']').toString();}
    private static float sf(float shape){float twice=2.0f*shape;return 1.0f<twice?1.0f:twice;}
    private static float inverse(float value,String relation,float shape){float n=value>ZERO?1.0f/value:INF;return relation.equals("prior")?sf(shape)*n:n;}
    private static void run(C c,boolean first){
        float[] mems=new float[c.types.length];float lesion=0.0f,bh=0.0f,relation=0.0f;
        float lesionEnergy=0.0f,bhEnergy=0.0f;
        for(int k=0;k<c.types.length;k++){
            String type=c.types[k],rel=c.relations[k];float shape=c.shapes[k];
            mems[k]=inverse(c.energies[k],rel,shape);
            if(type.equals("LESION")||type.equals("OPTIMIZED_LESION")){
                lesion=inverse(c.lesionE[k],rel,shape);relation=rel.equals("direct")?1.0f:sf(shape);
                lesionEnergy=c.lesionE[k];
                if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;relation*=n/c.maxVent;}
                else if(c.gmFlag&&c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;relation*=n/c.maxGM;}
                if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;relation*=(float)c.bstem/c.maxBstem;}
            }else if(type.equals("BLACKHOLE")){
                lesion=inverse(c.lesionE[k],rel,shape);bh=inverse(c.blackHoleE[k],rel,shape);relation=rel.equals("direct")?1.0f:sf(shape);
                lesionEnergy=c.lesionE[k];bhEnergy=c.blackHoleE[k];
                if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;relation*=n/c.maxVent;}
                else if(c.gmFlag&&c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;bh=n*bh/c.maxGM;relation*=n/c.maxGM;}
                if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;bh=c.bstem*bh/c.maxBstem;relation*=(float)c.bstem/c.maxBstem;}
            }
        }
        float den=0.0f;
        for(int k=0;k<c.types.length;k++){
            if(c.types[k].equals("LESION")||c.types[k].equals("OPTIMIZED_LESION"))den+=lesion;
            else if(c.types[k].equals("BLACKHOLE"))den+=lesion+bh;
            den+=mems[k];
        }
        float energy=0.0f;
        if(den>0.0f){for(int k=0;k<c.types.length;k++){
            mems[k]=mems[k]/den;
            if(c.types[k].equals("LESION")||c.types[k].equals("OPTIMIZED_LESION"))lesion=lesion/den;
            else if(c.types[k].equals("BLACKHOLE")){lesion=lesion/den;bh=bh/den;}
            energy+=c.energies[k]*mems[k]*mems[k];
        }}
        energy+=lesionEnergy*lesion*lesion+bhEnergy*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"types\":"+strings(c.types)+",\"relations\":"+strings(c.relations)
            +",\"class_energy_bits\":"+fbits(c.energies)+",\"shape_bits\":"+fbits(c.shapes)
            +",\"lesion_energy_bits\":"+fbits(c.lesionE)+",\"black_hole_energy_bits\":"+fbits(c.blackHoleE)
            +",\"ventricle_distance\":"+c.vent+",\"gray_matter_distance\":"+c.gm+",\"brainstem_distance\":"+c.bstem
            +",\"max_ventricle_distance\":"+c.maxVent+",\"max_gray_matter_distance\":"+c.maxGM+",\"max_brainstem_distance\":"+c.maxBstem
            +",\"gm_flag\":"+c.gmFlag+",\"spread_bits\":"+bits(c.spread)+",\"class_membership_bits\":"+fbits(mems)
            +",\"lesion_membership_bits\":"+bits(lesion)+",\"black_hole_membership_bits\":"+bits(bh)
            +",\"lesion_relation_bits\":"+bits(relation)+",\"denominator_bits\":"+bits(den)+",\"energy_bits\":"+bits(energy)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships selected formula composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);
        System.out.println("]}");
    }
}
