/** Formula excerpt for direct reciprocal-energy membership normalization. */
public final class ToadsDirectMembershipReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, CLASSES=4;
    private static final float ZERO=1e-10f, INF=1e10f;
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    private static float energy(int x,int y,int z,int c) {
        float value=(((13*x+7*y+3*z+11*c)%29)+1)/13.0f;
        if(x==2&&y==2&&z==2&&c==0) return 0.0f;
        if(x==3&&y==2&&z==2&&c==1) return ZERO;
        if(x==4&&y==3&&z==2&&c==2) return Math.nextAfter(ZERO,0.0f);
        if(x==2&&y==3&&z==2&&c==3) return Math.nextAfter(ZERO,1.0f);
        return value;
    }
    private static float objective(int x,int y,int z,int c) {
        return (((5*x+3*y+7*z+2*c)%19)+1)/17.0f;
    }
    private static float initial(int x,int y,int z,int c) {
        return ((5*x+7*y+11*z+3*c)%31)/30.0f;
    }
    private static float[] updated() {
        float[] output=new float[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            for(int c=0;c<CLASSES;c++) output[((x*NY+y)*NZ+z)*CLASSES+c]=initial(x,y,z,c);
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            float den=0.0f;
            for(int c=0;c<CLASSES;c++) {
                float num=energy(x,y,z,c);
                if(num>ZERO) num=1.0f/num;
                else num=INF;
                output[((x*NY+y)*NZ+z)*CLASSES+c]=num;
                den+=num;
            }
            if(den>0.0f) for(int c=0;c<CLASSES;c++) {
                int index=((x*NY+y)*NZ+z)*CLASSES+c;
                output[index]=output[index]/den;
            }
        }
        return output;
    }
    private static float[] contributions(float[] mems) {
        float[] output=new float[NX*NY*NZ*CLASSES];
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++)
            for(int c=0;c<CLASSES;c++) {
                int index=((x*NY+y)*NZ+z)*CLASSES+c;
                output[index]=objective(x,y,z,c)*mems[index]*mems[index];
            }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        float[] mems=updated();
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships direct reciprocal and normalization excerpt\",");
        System.out.print("\"source_sha256\":\""+SOURCE_SHA+"\",\"shape\":[7,6,5,4],\"zero\":1e-10,\"inf\":1e10,");
        System.out.print("\"updated\":");printArray(mems);
        System.out.print(",\"objective_contributions\":");printArray(contributions(mems));
        System.out.println('}');
    }
}
