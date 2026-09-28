/** Formula excerpt for source relation weighting plus reciprocal normalization. */
public final class ToadsShapeWeightedMembershipReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, CLASSES=4;
    private static final float PRIOR_FACTOR=0.75f, ZERO=1e-10f, INF=1e10f;
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String NUMERICS_SHA="13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";

    private static float energy(int x,int y,int z,int c) {
        float value=(((13*x+7*y+3*z+11*c)%29)+1)/13.0f;
        if(x==2&&y==2&&z==2&&c==0) return 0.0f;
        if(x==3&&y==2&&z==2&&c==1) return ZERO;
        if(x==4&&y==3&&z==2&&c==2) return Math.nextAfter(ZERO,0.0f);
        if(x==2&&y==3&&z==2&&c==3) return Math.nextAfter(ZERO,1.0f);
        return value;
    }
    private static float initial(int x,int y,int z,int c) {
        return ((5*x+7*y+11*z+3*c)%31)/30.0f;
    }
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
            for(int c=0;c<CLASSES;c++) output[((x*NY+y)*NZ+z)*CLASSES+c]=initial(x,y,z,c);
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            float den=0.0f;
            for(int c=0;c<CLASSES;c++) {
                boolean direct=label(x,y,z,0)==c||label(x,y,z,1)==c;
                boolean neighbor=label(x,y,z,2)==c||label(x,y,z,3)==c;
                float sh=shape(x,y,z,c), relationFactor;
                boolean priorOnly=!direct&&!neighbor;
                if(direct) relationFactor=1.0f;
                else if(neighbor) relationFactor=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));
                else if(sh==0.0f) relationFactor=0.0f;
                else relationFactor=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                float value;
                if(priorOnly&&sh==0.0f) value=0.0f;
                else {
                    float num=energy(x,y,z,c);
                    if(num>ZERO) num=1.0f/num;
                    else num=INF;
                    value=relationFactor*num;
                }
                output[((x*NY+y)*NZ+z)*CLASSES+c]=value;
                den+=value;
            }
            if(den>0.0f) for(int c=0;c<CLASSES;c++) {
                int index=((x*NY+y)*NZ+z)*CLASSES+c;
                output[index]=output[index]/den;
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
        System.out.print("{\"source\":\"MultipleToadSegmentation relation factors plus direct reciprocal normalization\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"numerics_source_sha256\":\""+NUMERICS_SHA+"\",");
        System.out.print("\"shape\":[7,6,5,4],\"prior_factor\":0.75,\"zero\":1e-10,\"inf\":1e10,\"updated\":");
        printArray(evaluate());
        System.out.println('}');
    }
}
