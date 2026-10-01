/** Pinned wrapper arguments and membership reference for updateShapeRegistration. */
public final class LesionToadsRegistrationSetupHandoffReferenceHarness {
    private static final class Atlas {
        float[][][][] memberships;
        int itMax, levels, subsample;
        float couplingFactor;
        String mode;
        void updateShapeRegistration(float[][][][] mems, String type, int iterations, int levelCount, float coupling) {
            memberships=mems;
            itMax=iterations;
            levels=levelCount;
            subsample=2;
            couplingFactor=coupling;
            mode=type;
        }
    }
    public static void main(String[] args) {
        float[][][][] membershipBuffer=new float[5][6][7][3];
        Atlas atlas=new Atlas();
        atlas.updateShapeRegistration(membershipBuffer,"rigid",1,1,1.0f);
        long sum=0;
        for(int x=0;x<5;x++)for(int y=0;y<6;y++)for(int z=0;z<7;z++)for(int k=0;k<3;k++)
            sum+=Float.floatToIntBits(atlas.memberships[x][y][z][k]);
        System.out.println("{\"source_method\":\"AlgorithmLesionToads to LesionToadAtlas.updateShapeRegistration\",\"shape\":[5,6,7,3],\"mode\":\""+atlas.mode+"\",\"iterations\":"+atlas.itMax+",\"levels\":"+atlas.levels+",\"subsample\":"+atlas.subsample+",\"coupling_factor_bits\":"+Integer.toUnsignedLong(Float.floatToIntBits(atlas.couplingFactor))+",\"same_membership_reference\":"+(atlas.memberships==membershipBuffer)+",\"membership_float_bits_sum\":"+sum+",\"call_order\":[\"segmentation.getMemberships\",\"atlas.updateShapeRegistration\"]}");
    }
}
