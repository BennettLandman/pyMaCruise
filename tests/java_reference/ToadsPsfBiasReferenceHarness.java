/** Source-formula excerpt of selected MultipleToadSegmentation PSF correction. */
public final class ToadsPsfBiasReferenceHarness {
    private static final int NX=9, NY=8, NZ=7, NC=2;
    private static final float RIDGE=0.37f;
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    private static float image(int x,int y,int z,int c) {
        return ((13*x+11*y+7*z+5*c)%113)/112.0f;
    }
    private static float field(int x,int y,int z,int c) {
        return 0.8f+((3*x+5*y+2*z+c)%21)/100.0f;
    }
    private static int label(int x,int y,int z) {
        if(x==2&&y==3&&z==3) return 2;
        if(x==6&&y>=1&&y<4&&z>=2&&z<5) return 0;
        return x>=4?1:0;
    }
    private static boolean isBoundary(int x,int y,int z) {
        int center=label(x,y,z);
        for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++)
            if(label(x+i,y+j,z+k)!=center) return true;
        return false;
    }
    private static float sample(int x,int y,int z,int c,boolean useField) {
        return useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
    }
    private static float[] evaluate(boolean useField) {
        float[] output=new float[NX*NY*NZ*NC];
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            if(!isBoundary(x,y,z)) continue;
            for(int c=0;c<NC;c++) {
                float value=-(27.0f/26.0f)*sample(x,y,z,c,useField);
                for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++)
                    value+=(1.0f/26.0f)*sample(x+i,y+j,z+k,c,useField);
                value*=RIDGE;
                output[((x*NY+y)*NZ+z)*NC+c]=value;
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
        System.out.print("{\"source\":\"MultipleToadSegmentation.isBoundary and computeExactMemberships PSF formula excerpt\",");
        System.out.print("\"source_sha256\":\""+SOURCE_SHA+"\",\"shape\":[9,8,7,2],\"ridge\":0.37,");
        System.out.print("\"without_field\":");printArray(evaluate(false));
        System.out.print(",\"with_field\":");printArray(evaluate(true));
        System.out.println('}');
    }
}
