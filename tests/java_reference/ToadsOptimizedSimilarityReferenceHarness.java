/** Formula excerpt of optimized-class pairwise similarity from computeExactMemberships. */
public final class ToadsOptimizedSimilarityReferenceHarness {
    private static final int NC=3, CLASSES=3;
    private static final float PRIOR_SCALE=0.05f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};
    // These are the no-field outputs of the committed P18am/P18an Java harnesses.
    private static final float[] VAR_IMG={0.41083512f,0.43782088f,0.15134391f};
    private static final float[][] STDDEV={{0.26655066f,0.25469658f,0.35447255f},
        {0.3371573f,0.2576227f,0.30383086f},{0.4673762f,0.29970098f,0.40097645f}};

    private static float[] evaluate() {
        float[] output=new float[CLASSES*CLASSES]; int index=0;
        for(int k=0;k<CLASSES;k++) for(int m=0;m<CLASSES;m++) {
            float similarity=0.0f;
            if(m!=k) {
                float centroidDist=0.0f;
                for(int c=0;c<NC;c++) {
                    centroidDist += VAR_IMG[c]*(CENTROIDS[c][k]-CENTROIDS[c][m])*(CENTROIDS[c][k]-CENTROIDS[c][m])
                        /(STDDEV[c][k]+STDDEV[c][m])/(STDDEV[c][k]+STDDEV[c][m]);
                }
                similarity=1.0f/(1.0f+(1.0f-PRIOR_SCALE)/PRIOR_SCALE*centroidDist);
            }
            output[index++]=similarity;
        }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships optimized centroid similarity formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"prior_scale\":0.05,");
        System.out.print("\"centroids\":[[0.64,0.62,0.81],[0.35,0.55,0.72],[0.18,0.7,0.88]],");
        System.out.print("\"image_variance_weights\":[0.41083512,0.43782088,0.15134391],");
        System.out.print("\"class_standard_deviations\":[[0.26655066,0.25469658,0.35447255],[0.3371573,0.2576227,0.30383086],[0.4673762,0.29970098,0.40097645]],");
        System.out.print("\"similarity\":");printArray(evaluate());System.out.println('}');
    }
}
