/** Source-formula excerpt of the six-face smoothing in computeExactMemberships. */
public final class ToadsSmoothingReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, CLASSES=4;
    private static final float SMOOTH=0.017f;
    private static final float[] PRIOR={0.5f,0.8f,1.1f,0.25f};
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";

    private static float membership(int x,int y,int z,int c) {
        return ((7*x+11*y+13*z+5*c)%47)/46.0f;
    }
    private static float[] evaluate(boolean energyTerm) {
        float[] output=new float[NX*NY*NZ*CLASSES];
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++)
            for(int k=0;k<CLASSES;k++) {
                float ngb=0.0f;
                for(int direction=0;direction<6;direction++) {
                    int xx=x, yy=y, zz=z;
                    if(direction==0) xx=x+1;
                    else if(direction==1) xx=x-1;
                    else if(direction==2) yy=y+1;
                    else if(direction==3) yy=y-1;
                    else if(direction==4) zz=z+1;
                    else zz=z-1;
                    for(int m=0;m<CLASSES;m++) if(m!=k) {
                        float value=membership(xx,yy,zz,m);
                        ngb+=value*value;
                    }
                }
                float value=energyTerm?0.5f*PRIOR[k]*SMOOTH*ngb:PRIOR[k]*SMOOTH*ngb;
                output[((x*NY+y)*NZ+z)*CLASSES+k]=value;
            }
        return output;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships six-face smoothing excerpt\",");
        System.out.print("\"source_sha256\":\""+SOURCE_SHA+"\",\"shape\":[7,6,5,4],\"smoothing\":0.017,");
        System.out.print("\"smoothing_priors\":[0.5,0.8,1.1,0.25],\"num\":");printArray(evaluate(false));
        System.out.print(",\"energy\":");printArray(evaluate(true));
        System.out.println('}');
    }
}
