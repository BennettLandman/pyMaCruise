/** Selected OUTMASK mixing through the shared voxel normalization loop. */
public final class LesionToadsOutmaskNormalizationReferenceHarness {
    private static final class C {
        final String[] types,relations;final float[] classValues,outmaskValues,shapes,energies;final float outlier;
        C(String[] t,String[] r,float[] c,float[] o,float[] s,float[] e,float l){types=t;relations=r;classValues=c;outmaskValues=o;shapes=s;energies=e;outlier=l;}
    }
    private static final C[] CASES={
        new C(new String[]{"OBJECT","OUTMASK"},new String[]{"direct","direct"},new float[]{0.3f,0.2f},new float[]{0.0f,0.4f},new float[]{0.5f,0.25f},new float[]{1.25f,2.5f},0.0f),
        new C(new String[]{"OBJECT","OUTMASK"},new String[]{"direct","prior"},new float[]{0.15f,0.35f},new float[]{0.0f,0.8f},new float[]{0.4f,0.3f},new float[]{0.5f,3.0f},1.0f),
        new C(new String[]{"OBJECT","OUTMASK","OUTMASK"},new String[]{"prior","direct","prior"},new float[]{0.125f,0.25f,0.375f},new float[]{0.0f,0.5f,0.75f},new float[]{0.2f,0.35f,0.45f},new float[]{0.75f,1.25f,2.0f},0.5f)
    };
    private static long bits(float x){return Float.floatToIntBits(x)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String strings(String[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append('"').append(a[i]).append('"');}return s.append(']').toString();}
    private static float sf(float shape){float twice=2.0f*shape;return 1.0f<twice?1.0f:twice;}
    private static void run(C c,boolean first){
        float[] mem=c.classValues.clone();float outlier=c.outlier;
        for(int k=0;k<c.types.length;k++)if(c.types[k].equals("OUTMASK")){
            float mixed=(1.0f-outlier)*mem[k]+outlier*c.outmaskValues[k];
            if(c.relations[k].equals("prior")){
                float factor=sf(c.shapes[k]);mixed=factor*mixed;outlier=factor*outlier;
            }
            mem[k]=mixed;outlier=outlier*c.outmaskValues[k];
        }
        float den=0;for(float value:mem)den+=value;
        float energy=0;
        if(den>0)for(int k=0;k<c.types.length;k++){
            mem[k]=mem[k]/den;
            if(c.types[k].equals("OUTMASK"))outlier=outlier/den;
            energy+=c.energies[k]*mem[k]*mem[k];
        }
        if(!first)System.out.print(',');
        System.out.print("{\"types\":"+strings(c.types)+",\"relations\":"+strings(c.relations)
            +",\"class_values_bits\":"+vec(c.classValues)+",\"outmask_values_bits\":"+vec(c.outmaskValues)
            +",\"shape_bits\":"+vec(c.shapes)+",\"energy_bits\":"+vec(c.energies)
            +",\"initial_outlier_bits\":"+bits(c.outlier)+",\"class_memberships_bits\":"+vec(mem)
            +",\"outlier_membership_bits\":"+bits(outlier)+",\"denominator_bits\":"+bits(den)+",\"local_energy_bits\":"+bits(energy)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships OUTMASK and normalization composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);System.out.println("]}");
    }
}
