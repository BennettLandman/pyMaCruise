/** Selected Java float32 reference for a composed OBJECT/OUTLIER membership sweep. */
public final class ToadsObjectMaskOutlierSweepReferenceHarness {
    private static final int NX=4, NY=4, NZ=4, CLASSES=2, OUTLIER_CLASS=1;
    private static final float SMOOTHING=0.13f, PRIOR_FACTOR=0.61f, ZERO=1e-10f, INF=1e10f;
    private static final float[] SMOOTH_PRIORS={0.7f,1.1f};
    private static float shape(int x,int y,int z,int k) { return 0.25f+0.01f*(x+2*y+3*z+k); }
    private static float initial(int x,int y,int z,int k) { return 0.12f+0.01f*(2*x+y+z+k); }
    private static float outlier(int x,int y,int z) { return 0.18f+0.01f*(x+2*y+z); }
    private static float energy(int x,int y,int z,int k) { return 0.8f+0.03f*(x+2*y+3*z+2*k); }
    private static float outlierEnergy(int x,int y,int z) { return 1.2f+0.02f*(2*x+y+z); }

    private static float[][][][] update(float[][][] outlierMap, float[] objectiveOut) {
        float[][][][] membership=new float[NX][NY][NZ][CLASSES];
        float[][][][] shapePrior=new float[NX][NY][NZ][CLASSES];
        float[][][][] data=new float[NX][NY][NZ][CLASSES];
        float[][][] outEnergy=new float[NX][NY][NZ];
        int[][][][] labels=new int[NX][NY][NZ][4];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            outlierMap[x][y][z]=outlier(x,y,z);outEnergy[x][y][z]=outlierEnergy(x,y,z);
            for(int k=0;k<CLASSES;k++) {
                membership[x][y][z][k]=initial(x,y,z,k);
                shapePrior[x][y][z][k]=shape(x,y,z,k);
                data[x][y][z][k]=energy(x,y,z,k);
            }
            for(int s=0;s<4;s++)labels[x][y][z][s]=-1;
        }
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        float total=0.0f;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++) {
            float denominator=0.0f;
            float[] objective={0.0f,0.0f};
            for(int k=0;k<CLASSES;k++) {
                float num=data[x][y][z][k],numO=outEnergy[x][y][z];
                float out=outlierMap[x][y][z];
                float obj=num;
                if(k==OUTLIER_CLASS)obj=(1.0f-out)*num+out*numO;
                float sum=0.0f;
                for(int[] d:offsets)for(int m=0;m<CLASSES;m++)if(m!=k) {
                    float value=membership[x+d[0]][y+d[1]][z+d[2]][m];
                    sum=sum+value*value;
                }
                float addition=SMOOTH_PRIORS[k]*SMOOTHING*sum;
                num=num+addition;
                if(k==OUTLIER_CLASS)numO=numO+addition;
                obj=obj+(0.5f*SMOOTH_PRIORS[k])*SMOOTHING*sum;
                boolean direct=false,neighbor=false;
                for(int s=0;s<2;s++)direct|=labels[x][y][z][s]==k;
                for(int s=2;s<4;s++)neighbor|=labels[x][y][z][s]==k;
                float factor=direct?1.0f:neighbor?Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*shapePrior[x][y][z][k],1.0f)):Math.min(1.0f,2.0f*PRIOR_FACTOR*shapePrior[x][y][z][k]);
                float reciprocal=num>ZERO?1.0f/num:INF;
                float value;
                if(k==OUTLIER_CLASS) {
                    float reciprocalO=numO>ZERO?1.0f/numO:INF;
                    value=factor*((1.0f-out)*reciprocal+out*reciprocalO);
                    outlierMap[x][y][z]=(factor*out)*reciprocalO;
                } else value=factor*reciprocal;
                membership[x][y][z][k]=value;
                denominator=denominator+value;
                objective[k]=obj;
            }
            if(denominator>0.0f)for(int k=0;k<CLASSES;k++) {
                membership[x][y][z][k]=membership[x][y][z][k]/denominator;
                if(k==OUTLIER_CLASS)outlierMap[x][y][z]=outlierMap[x][y][z]/denominator;
                total=total+objective[k]*membership[x][y][z][k]*membership[x][y][z][k];
            }
        }
        objectiveOut[0]=total;return membership;
    }
    private static void append(StringBuilder out,int[] values) { out.append('[');for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]);}out.append(']'); }
    private static int[] bits(float[][][][] a,int axis) { int[] b=new int[NX*NY*NZ*axis];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<axis;k++)b[p++]=Float.floatToIntBits(a[x][y][z][k]);return b; }
    private static int[] bits(float[][][] a) { int[] b=new int[NX*NY*NZ];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)b[p++]=Float.floatToIntBits(a[x][y][z]);return b; }
    private static int[] inputs(int kind) {
        int[] b=new int[NX*NY*NZ*((kind==1||kind==3)?1:CLASSES)];int p=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            if(kind==0)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(energy(x,y,z,k));
            else if(kind==1)b[p++]=Float.floatToIntBits(outlierEnergy(x,y,z));
            else if(kind==2)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(initial(x,y,z,k));
            else if(kind==3)b[p++]=Float.floatToIntBits(outlier(x,y,z));
            else if(kind==4)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(shape(x,y,z,k));
        }return b;
    }
    public static void main(String[] args) {
        float[][][] outlier=new float[NX][NY][NZ];float[] objective={0.0f};
        float[][][][] updated=update(outlier,objective);
        StringBuilder out=new StringBuilder("{\"shape\":[4,4,4],\"classes\":2,\"class_types\":[\"OBJECT\",\"OUTLIER\"]");
        out.append(",\"class_data_energy_bits\":");append(out,inputs(0));out.append(",\"outlier_data_energy_bits\":");append(out,inputs(1));
        out.append(",\"initial_memberships_bits\":");append(out,inputs(2));out.append(",\"initial_outlier_membership_bits\":");append(out,inputs(3));
        out.append(",\"shape_priors_bits\":");append(out,inputs(4));out.append(",\"closest_labels\":[");
        int[] labels=new int[NX*NY*NZ*4];java.util.Arrays.fill(labels,-1);append(out,labels);
        out.append("],\"smoothing_prior_bits\":");append(out,new int[]{Float.floatToIntBits(SMOOTH_PRIORS[0]),Float.floatToIntBits(SMOOTH_PRIORS[1])});
        out.append(",\"smoothing\":").append(SMOOTHING).append(",\"prior_factor\":").append(PRIOR_FACTOR);
        out.append(",\"updated_bits\":");append(out,bits(updated,CLASSES));out.append(",\"updated_outlier_bits\":");append(out,bits(outlier));
        out.append(",\"objective_bits\":").append(Float.floatToIntBits(objective[0])).append('}');System.out.println(out);
    }
}
