/** Formula excerpt of selected optimized-tissue computeExactImageVariance. */
public final class ToadsOptimizedVarianceReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=3;
    private static final float INF=1e10f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static float image(int x,int y,int z,int c) { return ((11*x+7*y+5*z+3*c)%101)/100.0f; }
    private static float field(int x,int y,int z,int c) { return 0.85f+((2*x+3*y+5*z+c)%16)/100.0f; }
    private static float membership(int x,int y,int z,int k) { return ((5*x+7*y+11*z+3*k)%31)/30.0f; }
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};

    private static float[][] evaluate(boolean useField) {
        float[] inverseVariance=new float[NC];
        float[] imageStddev=new float[NC];
        float sum=0.0f;
        for(int c=0;c<NC;c++) {
            float num=0.0f, den=0.0f;
            for(int k=0;k<CLASSES;k++) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                float mem=membership(x,y,z,k);
                float img=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float residual=img-CENTROIDS[c][k];
                float diff=residual*residual;
                num+=(mem*mem+1e-30f)*diff*diff;
                den+=(mem*mem+1e-30f);
            }
            if(num>0.0f) inverseVariance[c]=den/num;
            else inverseVariance[c]=INF;
            sum+=inverseVariance[c];
        }
        for(int c=0;c<NC;c++) {
            imageStddev[c]=(float)Math.sqrt(inverseVariance[c]);
            inverseVariance[c]=inverseVariance[c]/sum;
        }
        return new float[][]{inverseVariance,imageStddev};
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        float[][] noField=evaluate(false), withField=evaluate(true);
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactImageVariance optimized formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"shape\":[7,6,5,3],\"infinite_variance\":1e10,");
        System.out.print("\"without_field\":{\"normalized_inverse_variance\":");printArray(noField[0]);
        System.out.print(",\"image_stddev\":");printArray(noField[1]);
        System.out.print("},\"with_field\":{\"normalized_inverse_variance\":");printArray(withField[0]);
        System.out.print(",\"image_stddev\":");printArray(withField[1]);
        System.out.println("}}");
    }
}
