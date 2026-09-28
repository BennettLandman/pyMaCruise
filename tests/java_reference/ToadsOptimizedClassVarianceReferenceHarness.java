/** Formula excerpt of selected optimized-class computeExactVariances branch. */
public final class ToadsOptimizedClassVarianceReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=3;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};
    private static float image(int x,int y,int z,int c) { return ((11*x+7*y+5*z+3*c)%101)/100.0f; }
    private static float field(int x,int y,int z,int c) { return 0.85f+((2*x+3*y+5*z+c)%16)/100.0f; }
    private static float membership(int x,int y,int z,int k) { return ((5*x+7*y+11*z+3*k)%31)/30.0f; }
    private static int segmentation(int x,int y,int z) { return (x+2*y+z)%CLASSES; }

    private static float[] evaluate(boolean useField) {
        float[] output=new float[NC*CLASSES]; int index=0;
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f, den=0.0f;
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                if(segmentation(x,y,z)==k) {
                    float mem=membership(x,y,z,k);
                    float residual=(useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c))-CENTROIDS[c][k];
                    num+=(mem*mem+1e-30f)*residual*residual;
                    den+=(mem*mem+1e-30f);
                }
            }
            output[index++]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactVariances optimized formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"shape\":[7,6,5,3],");
        System.out.print("\"without_field\":");printArray(evaluate(false));
        System.out.print(",\"with_field\":");printArray(evaluate(true));
        System.out.println('}');
    }
}
