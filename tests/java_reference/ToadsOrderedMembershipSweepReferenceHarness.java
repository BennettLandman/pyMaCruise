/** Formula excerpt for one ordered TOADS direct/neighbor/prior membership sweep. */
public final class ToadsOrderedMembershipSweepReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, CLASSES=4;
    private static final float SMOOTH=0.017f, PRIOR_FACTOR=0.75f, ZERO=1e-10f, INF=1e10f;
    private static final float[] SMOOTH_PRIOR={0.5f,0.8f,1.1f,0.25f};
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String NUMERICS_SHA="13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";

    private static float dataEnergy(int x,int y,int z,int c) {
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
    private static float[] update(float[] total) {
        float[] mems=new float[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            for(int c=0;c<CLASSES;c++) mems[((x*NY+y)*NZ+z)*CLASSES+c]=initial(x,y,z,c);
        float dist=0.0f;
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            float den=0.0f;
            float[] energy=new float[CLASSES];
            float[] num=new float[CLASSES];
            boolean[] suppressed=new boolean[CLASSES];
            for(int c=0;c<CLASSES;c++) {
                num[c]=dataEnergy(x,y,z,c);
                energy[c]=num[c];
                if(SMOOTH>0.0f) {
                    float ngb=0.0f;
                    int[][] direction={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                    for(int[] d:direction) for(int m=0;m<CLASSES;m++) if(m!=c) {
                        int index=(((x+d[0])*NY+y+d[1])*NZ+z+d[2])*CLASSES+m;
                        ngb+=mems[index]*mems[index];
                    }
                    float smoothNum=SMOOTH_PRIOR[c]*SMOOTH*ngb;
                    num[c]+=smoothNum;
                    energy[c]+=0.5f*SMOOTH_PRIOR[c]*SMOOTH*ngb;
                }
            }
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
                if(priorOnly&&sh==0.0f) {
                    value=0.0f;
                    suppressed[c]=true;
                } else {
                    float reciprocal=num[c]>ZERO?1.0f/num[c]:INF;
                    value=relationFactor*reciprocal;
                }
                mems[((x*NY+y)*NZ+z)*CLASSES+c]=value;
                den+=value;
            }
            if(den>0.0f) for(int c=0;c<CLASSES;c++) {
                int index=((x*NY+y)*NZ+z)*CLASSES+c;
                mems[index]=mems[index]/den;
                dist+=energy[c]*mems[index]*mems[index];
            }
        }
        total[0]=dist;
        return mems;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        float[] total=new float[1], mems=update(total);
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships ordered smoothing and update excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"numerics_source_sha256\":\""+NUMERICS_SHA+"\",");
        System.out.print("\"shape\":[7,6,5,4],\"smoothing\":0.017,\"prior_factor\":0.75,\"smoothing_priors\":[0.5,0.8,1.1,0.25],");
        System.out.print("\"updated\":");printArray(mems);
        System.out.print(",\"objective\":"+Float.toString(total[0]));
        System.out.println('}');
    }
}
