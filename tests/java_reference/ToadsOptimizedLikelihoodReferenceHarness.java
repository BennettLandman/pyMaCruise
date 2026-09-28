import edu.jhmi.rad.medic.utilities.PowerTable;

/** Formula excerpt of selected MultipleToadSegmentation optimized tissue terms. */
public final class ToadsOptimizedLikelihoodReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=3;
    private static final String SEGMENTATION_SHA="739a770b866dcc76b8198d795a06f2ae90eecaa437bed50c2fe151bf72bf1fb1";
    private static final String POWER_TABLE_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final float[][] CENTROIDS={{0.22f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};
    private static final float[][] EXPONENTS={{0.80f,1.0f,1.0f},{1.0f,1.0f,1.0f},{1.0f,1.0f,1.0f}};
    private static final float[] MODWEIGHT={0.9f,1.2f,0.8f};
    private static final float[] VARIMG={0.75f,1.25f,0.9f};
    private static final float[][] VARWEIGHT={{1.0f,0.8f,1.1f},{0.25f,1.0f,0.7f},{1.1f,0.7f,0.9f}};

    private static float image(int x,int y,int z,int c) { return ((11*x+7*y+5*z+3*c)%101)/100.0f; }
    private static float field(int x,int y,int z,int c) { return 0.85f+((2*x+3*y+5*z+c)%16)/100.0f; }
    private static float bias(int x,int y,int z,int c) { return ((5*x+3*y+z+2*c)%17)/1000.0f; }
    private static String parameters(float[] values) {
        StringBuilder out=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.toString(values[i]));}
        return out.append(']').toString();
    }
    private static float[] evaluate(boolean useField) {
        PowerTable[][] tables=new PowerTable[NC][CLASSES];
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++)
            tables[c][k]=new PowerTable(0.0f,4.0f,0.000001f,EXPONENTS[c][k]);
        float[] output=new float[NX*NY*NZ*CLASSES]; int index=0;
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f;
            for(int c=0;c<NC;c++) {
                float corrected=(useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c))-bias(x,y,z,c);
                float residual=corrected-CENTROIDS[c][k];
                float distance=residual*residual;
                float weight=MODWEIGHT[c]*VARWEIGHT[c][k]*VARIMG[c];
                num+=weight*tables[c][k].lookup((double)distance,(double)EXPONENTS[c][k]);
            }
            output[index++]=num;
        }
        return output;
    }

    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation optimized CSF/GM/WM intensity-term excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",");
        System.out.print("\"power_table_source_sha256\":\""+POWER_TABLE_SHA+"\",");
        System.out.print("\"shape\":[7,6,5,3],\"centroids\":[[0.22,0.62,0.81],[0.35,0.55,0.72],[0.18,0.7,0.88]],");
        System.out.print("\"exponents\":[[0.8,1,1],[1,1,1],[1,1,1]],\"modality_weights\":[0.9,1.2,0.8],");
        System.out.print("\"image_variances\":[0.75,1.25,0.9],\"class_variance_weights\":[[1,0.8,1.1],[0.25,1,0.7],[1.1,0.7,0.9]],");
        System.out.print("\"without_field\":");printArray(evaluate(false));
        System.out.print(",\"with_field\":");printArray(evaluate(true));
        System.out.println('}');
    }
}
