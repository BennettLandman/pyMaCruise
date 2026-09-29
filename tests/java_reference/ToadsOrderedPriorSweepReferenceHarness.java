/** Formula excerpt composing data, smoothing, atlas prior and ordered update. */
public final class ToadsOrderedPriorSweepReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=4;
    private static final float SMOOTH=0.017f, PRIOR_FACTOR=0.75f, PRIOR_COEFFICIENT=0.15f, ZERO=1e-10f, INF=1e10f;
    private static final float[] SMOOTH_PRIOR={0.5f,0.8f,1.1f,0.25f};
    private static final float[] ATLAS_PRIOR={0.2f,0.35f,0.25f,0.45f};
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f,0.76f},{0.35f,0.55f,0.72f,0.42f},{0.18f,0.70f,0.88f,0.30f}};
    private static final float[][] STDDEV={{0.26655066f,0.25469658f,0.35447255f,0.31f},
        {0.3371573f,0.2576227f,0.30383086f,0.28f},{0.4673762f,0.29970098f,0.40097645f,0.36f}};
    private static final float[] VAR_IMG={0.41083512f,0.43782088f,0.15134391f};
    private static final float[][] SIMILARITY=similarities();

    private static float[][] similarities() {
        float[][] result=new float[CLASSES][CLASSES];
        float scaleFactor=(1.0f-0.05f)/0.05f;
        for(int k=0;k<CLASSES;k++) for(int m=0;m<CLASSES;m++) if(m!=k) {
            float distance=0.0f;
            for(int c=0;c<NC;c++) distance+=VAR_IMG[c]*(CENTROIDS[c][k]-CENTROIDS[c][m])*(CENTROIDS[c][k]-CENTROIDS[c][m])
                /(STDDEV[c][k]+STDDEV[c][m])/(STDDEV[c][k]+STDDEV[c][m]);
            result[k][m]=1.0f/(1.0f+scaleFactor*distance);
        }
        return result;
    }
    private static float dataEnergy(int x,int y,int z,int c) {
        float value=(((13*x+7*y+3*z+11*c)%29)+1)/13.0f;
        if(x==2&&y==2&&z==2&&c==0) return 0.0f;
        if(x==3&&y==2&&z==2&&c==1) return ZERO;
        if(x==4&&y==3&&z==2&&c==2) return Math.nextAfter(ZERO,0.0f);
        if(x==2&&y==3&&z==2&&c==3) return Math.nextAfter(ZERO,1.0f);
        return value;
    }
    private static float initial(int x,int y,int z,int c) { return ((5*x+7*y+11*z+3*c)%31)/30.0f; }
    private static float shape(int x,int y,int z,int c) { return ((x+2*y+3*z+4*c)%11)/10.0f; }
    private static int label(int x,int y,int z,int slot) {
        if(slot==0) return (x+y+z)%4;
        if(slot==1) return (x+y+z)%3==0?((x+y+z)%4+1)%4:-1;
        if(slot==2) return (2*x+y)%4;
        return (x+2*y+z+1)%4;
    }
    private static float[] update(float[] objectiveOut) {
        float[] mems=new float[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int c=0;c<CLASSES;c++)
            mems[((x*NY+y)*NZ+z)*CLASSES+c]=initial(x,y,z,c);
        float total=0.0f;
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            float denominator=0.0f;
            float[] energies=new float[CLASSES], num=new float[CLASSES];
            boolean[] suppressed=new boolean[CLASSES];
            for(int k=0;k<CLASSES;k++) {
                num[k]=dataEnergy(x,y,z,k); energies[k]=num[k];
                float neighborSum=0.0f;
                int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                for(int[] offset:offsets) for(int m=0;m<CLASSES;m++) if(m!=k) {
                    int index=(((x+offset[0])*NY+y+offset[1])*NZ+z+offset[2])*CLASSES+m;
                    neighborSum+=mems[index]*mems[index];
                }
                num[k]+=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                energies[k]+=0.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                float invprior=0.0f;
                for(int m=0;m<CLASSES;m++) if(m!=k)
                    invprior+=ATLAS_PRIOR[m]*SIMILARITY[k][m]*shape(x,y,z,m)*shape(x,y,z,m);
                num[k]+=PRIOR_COEFFICIENT*invprior;
                energies[k]+=0.5f*PRIOR_COEFFICIENT*invprior;
            }
            for(int k=0;k<CLASSES;k++) {
                boolean direct=label(x,y,z,0)==k||label(x,y,z,1)==k;
                boolean neighbor=label(x,y,z,2)==k||label(x,y,z,3)==k;
                float sh=shape(x,y,z,k), factor;
                boolean priorOnly=!direct&&!neighbor;
                if(direct) factor=1.0f;
                else if(neighbor) factor=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));
                else if(sh==0.0f) factor=0.0f;
                else factor=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                float value;
                if(priorOnly&&sh==0.0f) { value=0.0f; suppressed[k]=true; }
                else value=factor*(num[k]>ZERO?1.0f/num[k]:INF);
                mems[((x*NY+y)*NZ+z)*CLASSES+k]=value;
                denominator+=value;
            }
            if(denominator>0.0f) for(int k=0;k<CLASSES;k++) {
                int index=((x*NY+y)*NZ+z)*CLASSES+k;
                mems[index]=mems[index]/denominator;
                if(!suppressed[k]) total+=energies[k]*mems[index]*mems[index];
            }
        }
        objectiveOut[0]=total;
        return mems;
    }
    private static void printArray(float[] values) {
        System.out.print('[');
        for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}
        System.out.print(']');
    }
    public static void main(String[] args) {
        float[] objective=new float[1], updated=update(objective);
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships ordered smoothing and atlas-prior formula excerpt\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"shape\":[7,6,5,4],\"prior_coefficient\":0.15,");
        System.out.print("\"atlas_shape_priors\":[0.2,0.35,0.25,0.45],\"similarity\":[");
        for(int k=0;k<CLASSES;k++){if(k>0)System.out.print(',');printArray(SIMILARITY[k]);}
        System.out.print("],\"updated\":");printArray(updated);
        System.out.println(",\"objective\":"+Float.toString(objective[0])+"}");
    }
}
