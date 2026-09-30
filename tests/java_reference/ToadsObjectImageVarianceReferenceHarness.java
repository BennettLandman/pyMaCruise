/** Formula excerpt of ordinary-OBJECT computeExactImageVariance from the selected TOADS source. */
public final class ToadsObjectImageVarianceReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=4, NY=3, NZ=3, NC=2, CLASSES=3;
    private static final float[][] CENTROIDS={{0.5f,0.5f,0.5f},{0.31f,0.57f,0.76f}};

    private static float image(int x,int y,int z,int c) {
        if(c==0) return 0.5f;
        return (float)(0.12 + ((7*x+5*y+3*z)%19)*0.027);
    }
    private static float field(int x,int y,int z,int c) {
        if(c==0) return 1.0f;
        return (float)(0.8 + ((x+3*y+5*z)%7)*0.04);
    }
    private static float membership(int x,int y,int z,int k) {
        return (float)(((3*x+y+2*z+2*k)%7)/6.0);
    }
    private static int bits(float value) { return Float.floatToRawIntBits(value); }
    private static long u(float value) { return Integer.toUnsignedLong(bits(value)); }

    private static float[][] evaluate(boolean useField) {
        float[] numerator=new float[NC], denominator=new float[NC];
        float[] inverseVariance=new float[NC], stddev=new float[NC];
        float[] normalized=new float[NC];
        float sum=0.0f;
        for(int c=0;c<NC;c++) {
            float num=0.0f, den=0.0f;
            for(int k=0;k<CLASSES;k++) for(int x=0;x<NX;x++)
                for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                    float mem=membership(x,y,z,k);
                    float value=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                    float residual=value-CENTROIDS[c][k];
                    float diff=residual*residual;
                    num+=(mem*mem+1e-30f)*diff*diff;
                    den+=(mem*mem+1e-30f);
                }
            numerator[c]=num; denominator[c]=den;
            inverseVariance[c]=num>0.0f?den/num:1.0e10f;
            sum+=inverseVariance[c];
        }
        for(int c=0;c<NC;c++) {
            stddev[c]=(float)Math.sqrt(inverseVariance[c]);
            normalized[c]=inverseVariance[c]/sum;
        }
        return new float[][]{numerator,denominator,inverseVariance,{sum},stddev,normalized};
    }
    private static void print(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(u(values[i]));}
        System.out.print(']');
    }
    private static void printCase(String name, boolean useField) {
        float[][] r=evaluate(useField);
        System.out.print("\""+name+"\":{");
        System.out.print("\"numerator_bits\":");print(r[0]);
        System.out.print(",\"denominator_bits\":");print(r[1]);
        System.out.print(",\"inverse_variance_bits\":");print(r[2]);
        System.out.print(",\"inverse_variance_sum_bits\":");print(r[3]);
        System.out.print(",\"stddev_bits\":");print(r[4]);
        System.out.print(",\"normalized_var_img_bits\":");print(r[5]);
        System.out.print('}');
    }
    public static void main(String[] args) {
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"shape\":[4,3,3,2],\"class_count\":3,\"without_field\":{");
        float[][] noField=evaluate(false);
        System.out.print("\"numerator_bits\":");print(noField[0]);
        System.out.print(",\"denominator_bits\":");print(noField[1]);
        System.out.print(",\"inverse_variance_bits\":");print(noField[2]);
        System.out.print(",\"inverse_variance_sum_bits\":");print(noField[3]);
        System.out.print(",\"stddev_bits\":");print(noField[4]);
        System.out.print(",\"normalized_var_img_bits\":");print(noField[5]);
        System.out.print("},\"with_field\":{");
        float[][] withField=evaluate(true);
        System.out.print("\"numerator_bits\":");print(withField[0]);
        System.out.print(",\"denominator_bits\":");print(withField[1]);
        System.out.print(",\"inverse_variance_bits\":");print(withField[2]);
        System.out.print(",\"inverse_variance_sum_bits\":");print(withField[3]);
        System.out.print(",\"stddev_bits\":");print(withField[4]);
        System.out.print(",\"normalized_var_img_bits\":");print(withField[5]);
        System.out.println("}}");
    }
}
