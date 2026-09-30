/** Selected ordinary-class computeMemberships excerpt for GM ventricle scaling. */
public final class LesionGmVentricleMembershipReferenceHarness {
    private static final class C {
        final String name, relation;
        final String[] names;
        final float[] nums, shapes;
        final boolean hasLesion;
        final int ventDistance, maxVentDistance;
        C(String n, String r, String[] cn, float[] v, float[] sh, boolean lesion, int d, int max) {
            name=n; relation=r; names=cn; nums=v; shapes=sh; hasLesion=lesion;
            ventDistance=d; maxVentDistance=max;
        }
    }
    private static final C[] CASES = {
        new C("direct_exact_names_valid", "DIRECT",
            new String[]{"Cerebrum-GM","Other","GM","cerebrum-gm"},
            new float[]{0.3456789f,0.20000002f,0.1234567f,0.11111112f},
            new float[]{0.2f,0.4f,0.3f,0.1f}, true, 3, 8),
        new C("prior_exact_names_valid", "PRIOR",
            new String[]{"Cerebrum-GM","Other","GM","CerebrumGM"},
            new float[]{0.3456789f,0.20000002f,0.1234567f,0.11111112f},
            new float[]{0.2f,0.4f,0.3f,0.1f}, true, 3, 8),
        new C("distance_zero", "DIRECT",
            new String[]{"GM","Peer"}, new float[]{0.4f,0.6f}, new float[]{0.5f,0.5f}, true, 0, 8),
        new C("invalid_distance_falls_through", "DIRECT",
            new String[]{"GM","Peer"}, new float[]{0.4f,0.6f}, new float[]{0.5f,0.5f}, true, -1, 8),
        new C("zero_maximum_falls_through", "PRIOR",
            new String[]{"Cerebrum-GM","Peer"}, new float[]{0.4f,0.6f}, new float[]{0.25f,0.5f}, true, 2, 0),
        new C("no_lesion_falls_through", "PRIOR",
            new String[]{"GM","Peer"}, new float[]{0.4f,0.6f}, new float[]{0.25f,0.5f}, false, 2, 8),
        new C("near_miss_name", "PRIOR",
            new String[]{"CerebrumGM","cerebrum-gm","Gm","Peer"},
            new float[]{0.1f,0.2f,0.3f,0.4f}, new float[]{0.1f,0.2f,0.3f,0.4f}, true, 2, 8),
        new C("prior_float32_order", "PRIOR",
            new String[]{"GM","Peer","Cerebrum-GM"},
            new float[]{0.12345679f,0.2345679f,0.3456789f},
            new float[]{0.33333334f,0.7654321f,0.27182818f}, true, 5, 9),
        new C("nonpositive_denominator_is_not_normalized", "DIRECT",
            new String[]{"PeerA","PeerB"}, new float[]{-0.75f,0.25f},
            new float[]{0.5f,0.5f}, false, 2, 8)
    };
    private static long bits(float v) { return Float.floatToIntBits(v) & 0xffffffffL; }
    private static String q(String s) { return "\""+s+"\""; }
    private static String strings(String[] a) {
        StringBuilder b=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(q(a[i]));}
        return b.append(']').toString();
    }
    private static String ints(float[] a) {
        StringBuilder b=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}
        return b.append(']').toString();
    }
    private static void run(C c, boolean first) {
        float[] out=new float[c.nums.length]; float den=0.0f;
        boolean attenuate=c.hasLesion && c.maxVentDistance!=0 && c.ventDistance!=-1;
        for(int k=0;k<c.nums.length;k++){
            float value;
            if(c.relation.equals("DIRECT")) value=c.nums[k];
            else value=Math.min(1.0f,2.0f*c.shapes[k])*c.nums[k];
            if(attenuate && (c.names[k].equals("Cerebrum-GM")||c.names[k].equals("GM"))){
                if(c.relation.equals("DIRECT")) value=c.nums[k]*c.ventDistance/c.maxVentDistance;
                else value=Math.min(1.0f,2.0f*c.shapes[k])*c.nums[k]*c.ventDistance/c.maxVentDistance;
            }
            out[k]=value; den+=value;
        }
        if(den>0.0f) for(int k=0;k<out.length;k++) out[k]=out[k]/den;
        if(!first)System.out.print(',');
        System.out.print("{\"name\":"+q(c.name)+",\"relation\":"+q(c.relation)+",\"class_names\":"+strings(c.names)
            +",\"num_bits\":"+ints(c.nums)+",\"shape_bits\":"+ints(c.shapes)
            +",\"has_lesion\":"+c.hasLesion+",\"ventricle_distance\":"+c.ventDistance
            +",\"max_ventricle_distance\":"+c.maxVentDistance+",\"denominator_bits\":"+bits(den)
            +",\"normalized_membership_bits\":"+ints(out)+",\"normalized_sum_bits\":"+bits(sum(out))+"}");
    }
    private static float sum(float[] values) { float result=0.0f;for(float value:values)result+=value;return result; }
    public static void main(String[] args) {
        System.out.print("{\"source_method\":\"LesionToadSegmentation.computeMemberships ordinary GM class branch and normalization\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(CASES[i],i==0);
        System.out.println("]}");
    }
}
