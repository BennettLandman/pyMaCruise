/** Selected ordinary-object formula excerpt from MultipleToadSegmentation.computeExactVariances. */
public final class ToadsObjectVarianceReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=4, NY=3, NZ=3, NC=2, CLASSES=4;
    private static final float[][] CENTROID={{0.31f,0.73f,1.1f,1.37f},{0.43f,0.89f,1.2f,1.49f}};

    private static float image(int x,int y,int z,int c) {
        return 0.13f + ((7*x+5*y+3*z+11*c)%23)*0.017f;
    }
    private static float field(int x,int y,int z,int c) {
        return 0.71f + ((3*x+2*y+5*z+c)%9)*0.031f;
    }
    private static float membership(int x,int y,int z,int k) {
        if(k==1) return 0.0f;
        return ((3*x+y+2*z+k)%5)/4.0f;
    }
    private static int segmentation(int x,int y,int z) {
        if((x==0&&y==0&&z==0)||(x==NX-1&&y==NY-1&&z==NZ-1)) return 2;
        return (x+y+z)%2;
    }
    private static int bits(float value) { return Float.floatToRawIntBits(value); }

    private static float[] evaluate(boolean useField) {
        float[] out=new float[NC*CLASSES]; int index=0;
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f, den=0.0f;
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                if(segmentation(x,y,z)==k) {
                    float mem=membership(x,y,z,k);
                    float diff=(useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c))-CENTROID[c][k];
                    num+=(mem*mem+1e-30f)*diff*diff;
                    den+=(mem*mem+1e-30f);
                }
            }
            out[index++]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
        return out;
    }
    private static void printBits(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Integer.toUnsignedLong(bits(values[i])));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"shape\":[4,3,3,2],\"centroid_bits\":[");
        for(int c=0;c<NC;c++)for(int k=0;k<CLASSES;k++) {
            if(c>0||k>0)System.out.print(',');System.out.print(Integer.toUnsignedLong(bits(CENTROID[c][k])));
        }
        System.out.print("],\"without_field_bits\":");printBits(evaluate(false));
        System.out.print(",\"with_field_bits\":");printBits(evaluate(true));
        System.out.println('}');
    }
}
