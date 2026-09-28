/** Formula excerpt of direct/neighbor/prior gating from computeExactMemberships. */
public final class ToadsShapeRelationReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, CLASSES=4;
    private static final float PRIOR_FACTOR=0.75f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String NUMERICS_SHA="13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";

    private static float shape(int x,int y,int z,int c) {
        return ((x+2*y+3*z+4*c)%11)/10.0f;
    }
    private static int label(int x,int y,int z,int slot) {
        if(slot==0) return (x+y+z)%4;
        if(slot==1) return (x+y+z)%3==0?((x+y+z)%4+1)%4:-1;
        if(slot==2) return (2*x+y)%4;
        return (x+2*y+z+1)%4;
    }
    private static float[] evaluate() {
        float[] output=new float[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            for(int c=0;c<CLASSES;c++) {
                boolean direct=label(x,y,z,0)==c||label(x,y,z,1)==c;
                boolean neighbor=label(x,y,z,2)==c||label(x,y,z,3)==c;
                float sh=shape(x,y,z,c), value;
                if(direct) value=1.0f;
                else if(neighbor) value=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));
                else if(sh==0.0f) value=0.0f;
                else value=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                output[((x*NY+y)*NZ+z)*CLASSES+c]=value;
            }
        return output;
    }
    private static int[] suppressed() {
        int[] output=new int[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            for(int c=0;c<CLASSES;c++)
                output[((x*NY+y)*NZ+z)*CLASSES+c]=
                    (label(x,y,z,0)!=c&&label(x,y,z,1)!=c&&label(x,y,z,2)!=c&&label(x,y,z,3)!=c&&shape(x,y,z,c)==0.0f)?1:0;
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    private static void printArray(int[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(values[i]);}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation relation gate and Numerics shape factors\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"numerics_source_sha256\":\""+NUMERICS_SHA+"\",");
        System.out.print("\"shape\":[7,6,5,4],\"prior_factor\":0.75,\"factors\":");printArray(evaluate());
        System.out.print(",\"suppressed\":");printArray(suppressed());
        System.out.println('}');
    }
}
