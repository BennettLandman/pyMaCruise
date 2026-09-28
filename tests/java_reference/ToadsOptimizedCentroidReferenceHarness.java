import edu.jhmi.rad.medic.utilities.PowerTable;

/** Formula excerpt of selected optimized CSF/GM/WM damped centroid updates. */
public final class ToadsOptimizedCentroidReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=3;
    private static final float ZERO=1e-10f, SMOOTHNESS=0.1f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWER_TABLE_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};
    private static final float[][] EXPONENTS={{0.8f,1.0f,1.0f},{1.0f,1.0f,1.0f},{1.0f,1.0f,1.0f}};

    private static float image(int x,int y,int z,int c) { return ((11*x+7*y+5*z+3*c)%101)/100.0f; }
    private static float field(int x,int y,int z,int c) { return 0.85f+((2*x+3*y+5*z+c)%16)/100.0f; }
    private static float membership(int x,int y,int z,int k) { return ((5*x+7*y+11*z+3*k)%31)/30.0f; }
    private static float[] evaluate(boolean useField) {
        PowerTable[][] tables=new PowerTable[NC][CLASSES];
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++)
            tables[c][k]=new PowerTable(0.0f,4.0f,0.000001f,EXPONENTS[c][k]);
        float[] output=new float[NC*CLASSES]; int index=0;
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f, den=0.0f;
            for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
                float img=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float residual=img-CENTROIDS[c][k];
                float dist=residual*residual;
                float factor=1.0f;
                if(dist>ZERO) factor=(float)(tables[c][k].lookup((double)dist,(double)EXPONENTS[c][k])/dist);
                float mem=membership(x,y,z,k);
                num+=factor*mem*mem*img;
                den+=factor*mem*mem;
            }
            float updated=CENTROIDS[c][k];
            if(den>0.0f) updated=(1.0f-SMOOTHNESS)*num/den+SMOOTHNESS*updated;
            output[index++]=updated;
        }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactCentroids optimized damped formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"power_table_source_sha256\":\""+POWER_TABLE_SHA+"\",");
        System.out.print("\"shape\":[7,6,5,3],\"smoothness\":0.1,\"zero\":1e-10,");
        System.out.print("\"centroids\":[[0.64,0.62,0.81],[0.35,0.55,0.72],[0.18,0.7,0.88]],");
        System.out.print("\"exponents\":[[0.8,1,1],[1,1,1],[1,1,1]],\"without_field\":");printArray(evaluate(false));
        System.out.print(",\"with_field\":");printArray(evaluate(true));
        System.out.println('}');
    }
}
