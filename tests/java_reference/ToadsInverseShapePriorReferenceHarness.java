/** Formula excerpt of selected computeExactMemberships inverse shape prior. */
public final class ToadsInverseShapePriorReferenceHarness {
    private static final int NX=5, NY=6, NZ=4, CLASSES=3;
    private static final float PRIOR_COEFFICIENT=0.15f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float[] SHAPE_PRIORS={0.2f,0.35f,0.45f};
    private static final float[][] SIMILARITY={{0.0f,0.30548206f,0.16063479f},
        {0.30548206f,0.0f,0.3687615f},{0.16063479f,0.3687615f,0.0f}};
    private static float shape(int x,int y,int z,int k) { return ((3*x+5*y+7*z+11*k)%23)/22.0f; }
    private static float[] evaluate(boolean objective) {
        float[] output=new float[NX*NY*NZ*CLASSES]; int index=0;
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            for(int k=0;k<CLASSES;k++) {
                float value=0.0f;
                if(x>0&&x<NX-1&&y>0&&y<NY-1&&z>0&&z<NZ-1) {
                    float invprior=0.0f;
                    for(int m=0;m<CLASSES;m++) if(m!=k)
                        invprior += SHAPE_PRIORS[m]*SIMILARITY[k][m]*shape(x,y,z,m)*shape(x,y,z,m);
                    if(objective) value=0.5f*PRIOR_COEFFICIENT*invprior;
                    else value=PRIOR_COEFFICIENT*invprior;
                }
                output[index++]=value;
            }
        }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships inverse shape-prior formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"shape\":[5,6,4,3],\"prior_coefficient\":0.15,");
        System.out.print("\"shape_priors\":[0.2,0.35,0.45],\"similarity\":[[0.0,0.30548206,0.16063479],[0.30548206,0.0,0.3687615],[0.16063479,0.3687615,0.0]],");
        System.out.print("\"inverse_energy\":");printArray(evaluate(false));
        System.out.print(",\"objective_energy\":");printArray(evaluate(true));System.out.println('}');
    }
}
