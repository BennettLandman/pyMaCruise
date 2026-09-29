/** Selected lesionMembership/lesionRelation distance-factor excerpt. */
public final class LesionToadsDistanceModulationReferenceHarness {
    private static final float MEMBERSHIP=0.43f,RELATION=0.67f;
    private static final int[][] CASES={
        {0,4,2,9,5,5,1},
        {3,4,2,9,5,5,1},
        {-1,0,-1,9,5,5,1},
        {-1,3,-1,9,5,5,1},
        {-1,3,-1,9,5,5,0},
        {-1,-1,2,9,5,5,0},
        {-1,-1,0,9,5,5,0},
        {2,3,4,0,0,0,1},
        {0,3,-1,0,5,5,1},
        {-1,2,3,9,5,5,1},
        {-1,2,3,9,5,0,0},
        {4,0,1,7,0,4,1}
    };
    private static final float[] SPREADS={0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.375f};
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static void run(int i,boolean first) {
        int[] c=CASES[i];float spread=SPREADS[i];float membership=MEMBERSHIP,relation=RELATION;
        if(c[3]!=0&&c[0]!=-1) {
            float numerator=c[0]==0?spread:c[0];
            membership=numerator*membership/c[3];relation*=numerator/c[3];
        } else if(c[6]!=0&&c[4]!=0&&c[1]!=-1) {
            float numerator=c[1]==0?spread:c[1];
            membership=numerator*membership/c[4];relation*=numerator/c[4];
        }
        if(c[5]!=0&&c[2]!=-1) {
            membership=c[2]*membership/c[5];relation*=(float)c[2]/c[5];
        }
        if(!first)System.out.print(',');
        System.out.print("{\"ventricle_distance\":"+c[0]+",\"gray_matter_distance\":"+c[1]+",\"brainstem_distance\":"+c[2]+",\"max_ventricle_distance\":"+c[3]+",\"max_gray_matter_distance\":"+c[4]+",\"max_brainstem_distance\":"+c[5]+",\"gm_flag\":"+(c[6]!=0)+",\"spread_bits\":"+bits(spread)+",\"membership_bits\":"+bits(membership)+",\"relation_bits\":"+bits(relation)+"}");
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships lesion/lesionRelation distance modulation\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<CASES.length;i++)run(i,i==0);
        System.out.println("]}");
    }
}
