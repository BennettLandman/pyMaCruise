/** Selected OUTMASK class/outlier membership mixture excerpt. */
public final class LesionToadsOutmaskMembershipReferenceHarness {
    private static final float[][] CASES={
        {0.0f,0.2f,0.0f,0.0f},{0.25f,0.7f,0.0f,0.125f},{0.5f,0.4f,0.0f,0.25f},
        {1.0f,0.2f,0.0f,0.5f},{0.0f,0.9f,0.25f,0.0f},{0.25f,0.6f,0.25f,0.25f},
        {0.5f,0.3f,0.5f,0.5f},{1.0f,0.8f,0.5f,1.0f},{0.1f,0.75f,0.75f,0.4f},
        {0.375f,0.625f,1.0f,0.2f},{0.8765432f,0.1234567f,0.33333334f,0.73123455f},
        {0.2345679f,0.9876543f,0.6666667f,0.21987654f}
    };
    private static final String[] RELATIONS={"direct","direct","direct","direct","direct","direct",
        "prior","prior","prior","prior","prior","prior"};
    private static long bits(float x){return Float.floatToIntBits(x)&0xffffffffL;}
    private static void run(int i,boolean first){
        float[] c=CASES[i];float classValue=c[0],outmask=c[1],outlier=c[2],shape=c[3];
        float mixed=(1.0f-outlier)*classValue+outlier*outmask;
        if(RELATIONS[i].equals("prior")){
            float factor=Math.min(1.0f,2.0f*shape);
            mixed=factor*mixed;outlier=factor*outlier;
        }
        outlier=outlier*outmask;
        if(!first)System.out.print(',');
        System.out.print("{\"class_input_bits\":"+bits(c[0])+",\"outmask_input_bits\":"+bits(c[1])
            +",\"outlier_input_bits\":"+bits(c[2])+",\"shape_bits\":"+bits(shape)
            +",\"relation_type\":\""+RELATIONS[i]+"\",\"class_output_bits\":"+bits(mixed)
            +",\"outlier_output_bits\":"+bits(outlier)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships OUTMASK assignment\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(i,i==0);System.out.println("]}");
    }
}
