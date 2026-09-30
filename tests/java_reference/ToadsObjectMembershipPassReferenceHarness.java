/** Selected exact-membership composition for ordinary OBJECT classes. */
public final class ToadsObjectMembershipPassReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=4, NY=4, NZ=4, CHANNELS=2, CLASSES=3;
    private static final float SMOOTHING=0.17f, PRIOR_FACTOR=0.62f, PRIOR_COEFFICIENT=0.31f;
    private static final float ZERO=1e-10f, INF=1e10f;
    private static final float[] SMOOTH_PRIORS={0.6f,0.9f,1.2f};
    private static final float[] ATLAS_PRIORS={0.21f,0.34f,0.45f};
    private static final float[] MODALITY_WEIGHTS={0.75f,1.125f};
    private static final float[] IMAGE_VARIANCES={0.8f,1.25f};
    private static final float[][] CLASS_VARIANCE={{0.8f,1.1f,0.9f},{1.2f,0.7f,1.05f}};
    private static final float[][] SIMILARITY={{0.0f,0.4f,0.2f},{0.35f,0.0f,0.55f},{0.25f,0.45f,0.0f}};

    private static float image(int x,int y,int z,int c) { return (((7*x+5*y+3*z+4*c)%17)-8)/9.0f; }
    private static float field(int x,int y,int z,int c) { return 0.82f+((x+2*y+z+c)%5)*0.09f; }
    private static float bias(int x,int y,int z,int c) { return ((2*x+3*y+5*z+c)%9)/31.0f; }
    private static float initial(int x,int y,int z,int k) { return ((5*x+7*y+11*z+3*k)%31)/30.0f; }
    private static float shape(int x,int y,int z,int k) {
        int code=(x+2*y+3*z+4*k)%13;
        if(code==0) return 0.0f;
        if(code==1) return 0.07f;
        return code/13.0f;
    }
    private static int label(int x,int y,int z,int slot) {
        int code=x+2*y+3*z;
        if(slot==0) return code%CLASSES;
        if(slot==1) return code%4==0?(code+1)%CLASSES:-1;
        if(slot==2) return (2*x+y+z)%CLASSES;
        return (x+y+2*z+1)%CLASSES;
    }
    private static float[][][][] update(boolean useField,boolean useBias,float[] objectiveOut) {
        float[][][][] mems=new float[NX][NY][NZ][CLASSES];
        float[][][][] transformed=new float[NX][NY][NZ][CLASSES];
        int[][][][] closest=new int[NX][NY][NZ][4];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            for(int k=0;k<CLASSES;k++){mems[x][y][z][k]=initial(x,y,z,k);transformed[x][y][z][k]=shape(x,y,z,k);}
            for(int s=0;s<4;s++)closest[x][y][z][s]=label(x,y,z,s);
        }
        float[][][][] data=new float[NX][NY][NZ][CLASSES];
        float[][][][] priorNum=new float[NX][NY][NZ][CLASSES];
        float[][][][] priorObjective=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++) {
            float num=0.0f;
            for(int c=0;c<CHANNELS;c++) {
                float corrected=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                if(useBias)corrected=corrected-bias(x,y,z,c);
                float residual=corrected-((c*5+k*3-4)/10.0f);
                float weight=MODALITY_WEIGHTS[c]*CLASS_VARIANCE[c][k];weight=weight*IMAGE_VARIANCES[c];
                float contribution=weight*residual;contribution=contribution*residual;
                num=num+contribution;
            }
            data[x][y][z][k]=num;
        }
        float halfCoefficient=0.5f*PRIOR_COEFFICIENT;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<CLASSES;k++) {
            float inversePrior=0.0f;
            for(int m=0;m<CLASSES;m++)if(m!=k) {
                float term=ATLAS_PRIORS[m]*SIMILARITY[k][m];
                term=term*transformed[x][y][z][m];term=term*transformed[x][y][z][m];
                inversePrior=inversePrior+term;
            }
            priorNum[x][y][z][k]=PRIOR_COEFFICIENT*inversePrior;
            priorObjective[x][y][z][k]=halfCoefficient*inversePrior;
        }

        float total=0.0f;
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++) {
            float den=0.0f;float[] objective=new float[CLASSES];
            for(int k=0;k<CLASSES;k++) {
                float num=data[x][y][z][k],energy=num;
                if(SMOOTHING>0.0f) {
                    float ngb=0.0f;
                    for(int[] d:offsets)for(int m=0;m<CLASSES;m++)if(m!=k)
                        ngb=ngb+mems[x+d[0]][y+d[1]][z+d[2]][m]*mems[x+d[0]][y+d[1]][z+d[2]][m];
                    num=num+SMOOTH_PRIORS[k]*SMOOTHING*ngb;
                    energy=energy+0.5f*SMOOTH_PRIORS[k]*SMOOTHING*ngb;
                }
                num=num+priorNum[x][y][z][k];energy=energy+priorObjective[x][y][z][k];
                objective[k]=energy;
                boolean direct=closest[x][y][z][0]==k||closest[x][y][z][1]==k;
                boolean neighbor=closest[x][y][z][2]==k||closest[x][y][z][3]==k;
                float sh=transformed[x][y][z][k];boolean priorOnly=!direct&&!neighbor;
                float factor;
                if(direct)factor=1.0f;
                else if(neighbor)factor=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));
                else if(sh==0.0f)factor=0.0f;
                else factor=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                float value;
                if(priorOnly&&sh==0.0f)value=0.0f;
                else value=factor*(num>ZERO?1.0f/num:INF);
                mems[x][y][z][k]=value;den=den+value;
            }
            if(den>0.0f)for(int k=0;k<CLASSES;k++) {
                mems[x][y][z][k]=mems[x][y][z][k]/den;
                total=total+objective[k]*mems[x][y][z][k]*mems[x][y][z][k];
            }
        }
        objectiveOut[0]=total;return mems;
    }
    private static int[] floatBits(float[][][][] a,int finalAxis) {
        int[] out=new int[NX*NY*NZ*finalAxis];int p=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<finalAxis;k++)out[p++]=Float.floatToIntBits(a[x][y][z][k]);
        return out;
    }
    private static int[] floatBits2(float[][] a) {
        int[] out=new int[a.length*a[0].length];int p=0;for(float[] row:a)for(float v:row)out[p++]=Float.floatToIntBits(v);return out;
    }
    private static void append(StringBuilder out,int[] values) {
        out.append('[');for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]);}out.append(']');
    }
    private static float[][][][] images() {
        float[][][][] a=new float[NX][NY][NZ][CHANNELS];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=image(x,y,z,c);return a;
    }
    private static float[][] centroids() {
        float[][] a=new float[CHANNELS][CLASSES];for(int c=0;c<CHANNELS;c++)for(int k=0;k<CLASSES;k++)a[c][k]=((c*5+k*3-4)/10.0f);return a;
    }
    private static float[][][][] channelField(boolean useField) {
        if(!useField)return null;float[][][][] a=new float[NX][NY][NZ][CHANNELS];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=field(x,y,z,c);return a;
    }
    private static float[][][][] psfBias(boolean useBias) {
        if(!useBias)return null;float[][][][] a=new float[NX][NY][NZ][CHANNELS];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=bias(x,y,z,c);return a;
    }
    private static void appendVolume(StringBuilder out,float[][][][] a,int axis) {
        if(a==null){out.append("null");return;}int[] bits=new int[NX*NY*NZ*axis];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<axis;c++)bits[p++]=Float.floatToIntBits(a[x][y][z][c]);append(out,bits);
    }
    private static void emit(StringBuilder out,String name,boolean useField,boolean useBias) {
        float[] objective={0.0f};float[][][][] updated=update(useField,useBias,objective);
        if(out.charAt(out.length()-1)!='[')out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"intensity_field\":");appendVolume(out,channelField(useField),CHANNELS);
        out.append(",\"psf_bias\":");appendVolume(out,psfBias(useBias),CHANNELS);
        out.append(",\"updated_bits\":");append(out,floatBits(updated,CLASSES));
        out.append(",\"objective_bits\":").append(Float.floatToIntBits(objective[0])).append('}');
    }
    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[4,4,4,2],\"classes\":3");
        out.append(",\"images_bits\":");appendVolume(out,images(),CHANNELS);out.append(",\"centroids_bits\":");append(out,floatBits2(centroids()));
        out.append(",\"modality_weights_bits\":");append(out,new int[]{Float.floatToIntBits(MODALITY_WEIGHTS[0]),Float.floatToIntBits(MODALITY_WEIGHTS[1])});
        out.append(",\"class_variance_bits\":");append(out,floatBits2(CLASS_VARIANCE));
        out.append(",\"image_variances_bits\":");append(out,new int[]{Float.floatToIntBits(IMAGE_VARIANCES[0]),Float.floatToIntBits(IMAGE_VARIANCES[1])});
        out.append(",\"initial_memberships_bits\":[");
        int[] initialBits=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)initialBits[p++]=Float.floatToIntBits(initial(x,y,z,k));append(out,initialBits);
        out.append("],\"transformed_shapes_bits\":[");int[] shapeBits=new int[NX*NY*NZ*CLASSES];p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)shapeBits[p++]=Float.floatToIntBits(shape(x,y,z,k));append(out,shapeBits);
        out.append("],\"closest_labels\":[");int[] labelValues=new int[NX*NY*NZ*4];p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int s=0;s<4;s++)labelValues[p++]=label(x,y,z,s);append(out,labelValues);
        out.append("],\"similarity_bits\":");append(out,floatBits2(SIMILARITY));
        out.append(",\"smoothing_prior_bits\":");append(out,new int[]{Float.floatToIntBits(SMOOTH_PRIORS[0]),Float.floatToIntBits(SMOOTH_PRIORS[1]),Float.floatToIntBits(SMOOTH_PRIORS[2])});
        out.append(",\"atlas_prior_bits\":");append(out,new int[]{Float.floatToIntBits(ATLAS_PRIORS[0]),Float.floatToIntBits(ATLAS_PRIORS[1]),Float.floatToIntBits(ATLAS_PRIORS[2])});
        out.append(",\"smoothing\":").append(SMOOTHING).append(",\"prior_factor\":").append(PRIOR_FACTOR)
            .append(",\"prior_coefficient\":").append(PRIOR_COEFFICIENT);
        out.append(",\"cases\":[");emit(out,"raw_no_bias",false,false);emit(out,"field_with_psf_bias",true,true);
        System.out.println(out.append("]}"));
    }
}
