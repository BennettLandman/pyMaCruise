/** Selected final one-voxel computeMemberships normalization excerpt. */
public final class LesionToadsVoxelMembershipNormalizationReferenceHarness {
    private static final class C {
        final String[] types; final float[] memberships,energies;
        final float lesion,blackHole,outlier,lesionEnergy,blackHoleEnergy;
        C(String[] t,float[] m,float[] e,float l,float b,float o,float le,float be){
            types=t;memberships=m;energies=e;lesion=l;blackHole=b;outlier=o;
            lesionEnergy=le;blackHoleEnergy=be;
        }
    }
    private static final C[] CASES={
        new C(new String[]{"OBJECT","OBJECT"},new float[]{0.3f,0.7f},new float[]{1.25f,2.0f},0.4f,0.2f,0.1f,3.0f,4.0f),
        new C(new String[]{"OBJECT","LESION"},new float[]{0.3f,0.2f},new float[]{1.5f,2.5f},0.4f,0.0f,0.2f,3.25f,0.0f),
        new C(new String[]{"OBJECT","BLACKHOLE"},new float[]{0.3f,0.2f},new float[]{1.5f,2.5f},0.4f,0.1f,0.2f,3.25f,4.5f),
        new C(new String[]{"OBJECT","OUTLIER","OUTMASK"},new float[]{0.2f,0.3f,0.1f},new float[]{0.75f,1.25f,2.5f},0.4f,0.2f,0.6f,3.0f,4.0f),
        new C(new String[]{"OBJECT","LESION"},new float[]{-0.5f,-0.5f},new float[]{1.0f,2.0f},-0.25f,0.0f,0.4f,3.0f,0.0f),
        new C(new String[]{"LESION","LESION"},new float[]{0.1f,0.2f},new float[]{0.75f,1.5f},0.3f,0.0f,0.2f,2.0f,0.0f),
        new C(new String[]{"OPTIMIZED_LESION","OBJECT"},new float[]{0.125f,0.375f},new float[]{4.0f,0.125f},0.25f,0.0f,0.1f,0.5f,0.0f),
        new C(new String[]{"MASK","OBJECT","BLACKHOLE"},new float[]{0.10000001f,0.20000002f,0.30000004f},new float[]{0.33333334f,1.2345678f,9.876543f},0.07123457f,0.01987654f,0.49999997f,2.345678f,7.654321f)
    };
    private static long bits(float v){return Float.floatToIntBits(v)&0xffffffffL;}
    private static String ints(float[] a){
        StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();
    }
    private static String strings(String[] a){
        StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('"').append(a[i]).append('"');}return b.append(']').toString();
    }
    private static void run(C c,boolean first){
        float[] m=c.memberships.clone();float lesion=c.lesion,bh=c.blackHole,outlier=c.outlier;
        float den=0.0f;
        for(int k=0;k<c.types.length;k++){
            if(c.types[k].equals("LESION")||c.types[k].equals("OPTIMIZED_LESION")) den+=lesion;
            else if(c.types[k].equals("BLACKHOLE")) den+=lesion+bh;
            den+=m[k];
        }
        float dist=0.0f;
        if(den>0.0f){
            for(int k=0;k<c.types.length;k++){
                m[k]=m[k]/den;
                if(c.types[k].equals("LESION")||c.types[k].equals("OPTIMIZED_LESION")) lesion=lesion/den;
                else if(c.types[k].equals("BLACKHOLE")){lesion=lesion/den;bh=bh/den;}
                else if(c.types[k].equals("OUTLIER")||c.types[k].equals("OUTMASK")) outlier=outlier/den;
                dist+=c.energies[k]*m[k]*m[k];
            }
        }
        dist+=c.lesionEnergy*lesion*lesion+c.blackHoleEnergy*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"types\":"+strings(c.types)+",\"class_membership_bits\":"+ints(c.memberships)
            +",\"class_energy_bits\":"+ints(c.energies)+",\"lesion_membership_bits\":"+bits(c.lesion)
            +",\"black_hole_membership_bits\":"+bits(c.blackHole)+",\"outlier_membership_bits\":"+bits(c.outlier)
            +",\"lesion_energy_bits\":"+bits(c.lesionEnergy)+",\"black_hole_energy_bits\":"+bits(c.blackHoleEnergy)
            +",\"denominator_bits\":"+bits(den)+",\"normalized_class_membership_bits\":"+ints(m)
            +",\"normalized_lesion_membership_bits\":"+bits(lesion)+",\"normalized_black_hole_membership_bits\":"+bits(bh)
            +",\"normalized_outlier_membership_bits\":"+bits(outlier)+",\"energy_bits\":"+bits(dist)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships denominator/normalization excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);
        System.out.println("]}");
    }
}
