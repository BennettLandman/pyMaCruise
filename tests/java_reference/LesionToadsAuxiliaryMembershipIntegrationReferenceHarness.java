/** Selected source composition of lesion/black-hole and OUTMASK auxiliaries. */
public final class LesionToadsAuxiliaryMembershipIntegrationReferenceHarness {
    private static final class C {
        final String auxiliary,relation;final float[] memberships,energies;final float lesion,blackHole,lesionRelation,outlier,lesionEnergy,bhEnergy;
        final float outmaskValue,outmaskShape;final String outmaskRelation;final int vent,gm,bstem,maxVent,maxGM,maxBstem;final float spread;
        C(String a,String r,float[] m,float[] e,float l,float b,float lr,float o,float le,float be,float ov,float os,String or,
          int v,int g,int bs,int mv,int mg,int mb,float sp){
            auxiliary=a;relation=r;memberships=m;energies=e;lesion=l;blackHole=b;lesionRelation=lr;outlier=o;lesionEnergy=le;bhEnergy=be;
            outmaskValue=ov;outmaskShape=os;outmaskRelation=or;vent=v;gm=g;bstem=bs;maxVent=mv;maxGM=mg;maxBstem=mb;spread=sp;
        }
    }
    private static final C[] CASES={
        new C("LESION","prior",new float[]{0.3f,0.2f,0.1f},new float[]{0.5f,1.25f,0.75f},0.25f,0.0f,0.6f,0.2f,2.0f,0.0f,0.4f,0.25f,"prior",-1,0,2,9,5,5,0.375f),
        new C("BLACKHOLE","direct",new float[]{0.35f,0.15f,0.2f},new float[]{0.6f,1.5f,0.8f},0.3f,0.4f,1.0f,0.75f,1.25f,2.5f,0.55f,0.4f,"direct",2,0,-1,4,5,5,0.5f)
    };
    private static long bits(float x){return Float.floatToIntBits(x)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String strings(String[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append('"').append(a[i]).append('"');}return s.append(']').toString();}
    private static float sf(float shape){float twice=2.0f*shape;return 1.0f<twice?1.0f:twice;}
    private static void run(C c,boolean first){
        String[] types={"OBJECT",c.auxiliary,"OUTMASK"};float[] mem=c.memberships.clone();
        float lesion=c.lesion,bh=c.blackHole,relation=c.lesionRelation,outlier=c.outlier;
        if(c.auxiliary.equals("LESION")){
            if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;if(c.vent==0)relation*=n/c.maxVent;else relation*=c.vent/c.maxVent;}
            else if(c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;if(c.gm==0)relation*=n/c.maxGM;else relation*=c.gm/c.maxGM;}
            if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;relation*=c.bstem/c.maxBstem;}
        }else{
            if(c.maxVent!=0&&c.vent!=-1){float n=c.vent==0?c.spread:c.vent;lesion=n*lesion/c.maxVent;if(c.vent==0)relation*=n/c.maxVent;else relation*=c.vent/c.maxVent;}
            else if(c.maxGM!=0&&c.gm!=-1){float n=c.gm==0?c.spread:c.gm;lesion=n*lesion/c.maxGM;bh=n*bh/c.maxGM;if(c.gm==0)relation*=n/c.maxGM;else relation*=c.gm/c.maxGM;}
            if(c.maxBstem!=0&&c.bstem!=-1){lesion=c.bstem*lesion/c.maxBstem;bh=c.bstem*bh/c.maxBstem;relation*=c.bstem/c.maxBstem;}
        }
        float classOutmask=(1.0f-outlier)*mem[2]+outlier*c.outmaskValue;
        if(c.outmaskRelation.equals("prior")){float factor=sf(c.outmaskShape);classOutmask=factor*classOutmask;outlier=factor*outlier;}
        mem[2]=classOutmask;outlier=outlier*c.outmaskValue;
        float den=0;for(int k=0;k<types.length;k++){if(types[k].equals("LESION"))den+=lesion;else if(types[k].equals("BLACKHOLE"))den+=lesion+bh;den+=mem[k];}
        float energy=0;if(den>0)for(int k=0;k<types.length;k++){
            mem[k]/=den;
            if(types[k].equals("LESION"))lesion/=den;
            else if(types[k].equals("BLACKHOLE")){lesion/=den;bh/=den;}
            if(types[k].equals("OUTMASK"))outlier/=den;
            energy+=c.energies[k]*mem[k]*mem[k];
        }
        energy+=c.lesionEnergy*lesion*lesion+c.bhEnergy*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"auxiliary_type\":\""+c.auxiliary+"\",\"distance_relation\":\""+c.relation+"\",\"class_membership_bits\":"+vec(c.memberships)
            +",\"class_energy_bits\":"+vec(c.energies)+",\"lesion_bits\":"+bits(c.lesion)+",\"black_hole_bits\":"+bits(c.blackHole)
            +",\"lesion_relation_bits\":"+bits(c.lesionRelation)+",\"outlier_bits\":"+bits(c.outlier)
            +",\"lesion_energy_bits\":"+bits(c.lesionEnergy)+",\"black_hole_energy_bits\":"+bits(c.bhEnergy)
            +",\"outmask_value_bits\":"+bits(c.outmaskValue)+",\"outmask_shape_bits\":"+bits(c.outmaskShape)
            +",\"outmask_relation\":\""+c.outmaskRelation+"\",\"ventricle_distance\":"+c.vent+",\"gray_matter_distance\":"+c.gm
            +",\"brainstem_distance\":"+c.bstem+",\"max_ventricle_distance\":"+c.maxVent+",\"max_gray_matter_distance\":"+c.maxGM
            +",\"max_brainstem_distance\":"+c.maxBstem+",\"spread_bits\":"+bits(c.spread)
            +",\"class_memberships_output_bits\":"+vec(mem)+",\"lesion_output_bits\":"+bits(lesion)+",\"black_hole_output_bits\":"+bits(bh)
            +",\"lesion_relation_output_bits\":"+bits(relation)+",\"outlier_output_bits\":"+bits(outlier)
            +",\"denominator_bits\":"+bits(den)+",\"local_energy_bits\":"+bits(energy)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships auxiliary/OUTMASK composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);System.out.println("]}");
    }
}
