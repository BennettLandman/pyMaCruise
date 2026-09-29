import edu.jhmi.rad.medic.utilities.PowerTable;

/** Source-formula harness for one optimized-tissue membership pass. */
public final class ToadsOptimizedMembershipPassReferenceHarness {
    private static final int NX=7, NY=6, NZ=5, NC=3, CLASSES=3;
    private static final float ZERO=1e-10f, INF=1e10f, SMOOTH=0.017f;
    private static final float PRIOR_FACTOR=0.75f, PRIOR_COEFFICIENT=0.15f, PRIOR_SCALE=0.05f;
    private static final float[] SMOOTH_PRIOR={0.5f,0.8f,1.1f};
    private static final float[] ATLAS_PRIOR={0.2f,0.35f,0.45f};
    private static final float[] MODALITY_WEIGHT={1.0f,0.75f,0.5f};
    private static final float[][] CENTROIDS={{0.64f,0.62f,0.81f},{0.35f,0.55f,0.72f},{0.18f,0.70f,0.88f}};
    private static final float[][] EXPONENTS={{0.8f,1.0f,1.0f},{1.0f,1.0f,1.0f},{1.0f,1.0f,1.0f}};
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static float image(int x,int y,int z,int c) { return ((11*x+7*y+5*z+3*c)%101)/100.0f; }
    private static float field(int x,int y,int z,int c) { return 0.85f+((2*x+3*y+5*z+c)%16)/100.0f; }
    private static float membership(int x,int y,int z,int k) { return ((5*x+7*y+11*z+3*k)%31)/30.0f; }
    private static float shape(int x,int y,int z,int k) { return ((x+2*y+3*z+4*k)%11)/10.0f; }
    private static int segmentation(int x,int y,int z) { return (x+2*y+z)%CLASSES; }
    private static int closest(int x,int y,int z,int slot) {
        int label;
        if(slot==0) label=(x+y+z)%4;
        else if(slot==1) label=(x+y+z)%3==0?((x+y+z)%4+1)%4:-1;
        else if(slot==2) label=(2*x+y)%4;
        else label=(x+2*y+z+1)%4;
        return label<CLASSES?label:-1;
    }
    private static float[][] imageVariance(boolean useField) {
        float[] inverse=new float[NC], stddev=new float[NC]; float sum=0.0f;
        for(int c=0;c<NC;c++) {
            float num=0.0f,den=0.0f;
            for(int k=0;k<CLASSES;k++) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                float mem=membership(x,y,z,k);
                float value=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float residual=value-CENTROIDS[c][k], diff=residual*residual;
                num+=(mem*mem+1e-30f)*diff*diff;
                den+=(mem*mem+1e-30f);
            }
            inverse[c]=num>0.0f?den/num:INF; sum+=inverse[c];
        }
        for(int c=0;c<NC;c++){stddev[c]=(float)Math.sqrt(inverse[c]);inverse[c]=inverse[c]/sum;}
        return new float[][]{inverse,stddev};
    }
    private static float[][] classStddev(boolean useField) {
        float[][] output=new float[NC][CLASSES];
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f,den=0.0f;
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) if(segmentation(x,y,z)==k) {
                float mem=membership(x,y,z,k);
                float value=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float diff=value-CENTROIDS[c][k];
                num+=(mem*mem+1e-30f)*diff*diff;
                den+=(mem*mem+1e-30f);
            }
            output[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
        return output;
    }
    private static float[][] similarities(float[] weights,float[][] deviations) {
        float[][] output=new float[CLASSES][CLASSES];
        for(int k=0;k<CLASSES;k++) for(int m=0;m<CLASSES;m++) if(m!=k) {
            float distance=0.0f;
            for(int c=0;c<NC;c++) distance+=weights[c]*(CENTROIDS[c][k]-CENTROIDS[c][m])*(CENTROIDS[c][k]-CENTROIDS[c][m])
                /(deviations[c][k]+deviations[c][m])/(deviations[c][k]+deviations[c][m]);
            output[k][m]=1.0f/(1.0f+(1.0f-PRIOR_SCALE)/PRIOR_SCALE*distance);
        }
        return output;
    }
    private static float[][][][] bias(boolean useField,float ridge) {
        float[][][][] output=new float[NX][NY][NZ][NC];
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            boolean boundary=false;
            for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++)
                if(segmentation(x+i,y+j,z+l)!=segmentation(x,y,z)) boundary=true;
            if(boundary) for(int c=0;c<NC;c++) {
                float center=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float value=-(27.0f/26.0f)*center;
                for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++) {
                    float sample=useField?field(x+i,y+j,z+l,c)*image(x+i,y+j,z+l,c):image(x+i,y+j,z+l,c);
                    value+=(1.0f/26.0f)*sample;
                }
                output[x][y][z][c]=value*ridge;
            }
        }
        return output;
    }
    private static float[][][][] dataEnergy(boolean useField,float[] varImg,float[][][][] psfBias) {
        PowerTable[][] tables=new PowerTable[NC][CLASSES];
        for(int c=0;c<NC;c++) for(int k=0;k<CLASSES;k++) tables[c][k]=new PowerTable(0.0f,4.0f,0.000001f,EXPONENTS[c][k]);
        float[][][][] output=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int k=0;k<CLASSES;k++) {
            float num=0.0f;
            for(int c=0;c<NC;c++) {
                float value=useField?field(x,y,z,c)*image(x,y,z,c):image(x,y,z,c);
                float residual=value-psfBias[x][y][z][c]-CENTROIDS[c][k];
                float distance=residual*residual;
                float weight=MODALITY_WEIGHT[c]*1.0f*varImg[c];
                num+=weight*(float)tables[c][k].lookup((double)distance,(double)EXPONENTS[c][k]);
            }
            output[x][y][z][k]=num;
        }
        return output;
    }
    private static float[] update(float[][][][] data,float[][][][] bias,float[][] similarity,float[] imageWeights,float[][] deviations) {
        float[] mems=new float[NX*NY*NZ*CLASSES];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) for(int k=0;k<CLASSES;k++)
            mems[((x*NY+y)*NZ+z)*CLASSES+k]=membership(x,y,z,k);
        float total=0.0f;
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            float den=0.0f; float[] num=new float[CLASSES], energy=new float[CLASSES]; boolean[] suppressed=new boolean[CLASSES];
            for(int k=0;k<CLASSES;k++) {
                num[k]=data[x][y][z][k]; energy[k]=num[k];
                float ngb=0.0f; int[][] dirs={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                for(int[] d:dirs) for(int m=0;m<CLASSES;m++) if(m!=k) {
                    float v=mems[(((x+d[0])*NY+y+d[1])*NZ+z+d[2])*CLASSES+m]; ngb+=v*v;
                }
                num[k]+=SMOOTH_PRIOR[k]*SMOOTH*ngb;
                energy[k]+=0.5f*SMOOTH_PRIOR[k]*SMOOTH*ngb;
                float inv=0.0f;
                for(int m=0;m<CLASSES;m++) if(m!=k) inv+=ATLAS_PRIOR[m]*similarity[k][m]*shape(x,y,z,m)*shape(x,y,z,m);
                num[k]+=PRIOR_COEFFICIENT*inv; energy[k]+=0.5f*PRIOR_COEFFICIENT*inv;
            }
            for(int k=0;k<CLASSES;k++) {
                boolean direct=closest(x,y,z,0)==k||closest(x,y,z,1)==k;
                boolean neighbor=closest(x,y,z,2)==k||closest(x,y,z,3)==k;
                float sh=shape(x,y,z,k), factor; boolean priorOnly=!direct&&!neighbor;
                if(direct) factor=1.0f;
                else if(neighbor) factor=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));
                else if(sh==0.0f) factor=0.0f;
                else factor=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                float v;
                if(priorOnly&&sh==0.0f){v=0.0f;suppressed[k]=true;}
                else v=factor*(num[k]>ZERO?1.0f/num[k]:INF);
                mems[((x*NY+y)*NZ+z)*CLASSES+k]=v; den+=v;
            }
            if(den>0.0f) for(int k=0;k<CLASSES;k++) {
                int idx=((x*NY+y)*NZ+z)*CLASSES+k; mems[idx]=mems[idx]/den;
                if(!suppressed[k]) total+=energy[k]*mems[idx]*mems[idx];
            }
        }
        return new float[]{total};
    }
    private static void printArray(float[] values) {
        System.out.print('[');for(int i=0;i<values.length;i++){if(i>0)System.out.print(',');System.out.print(Float.toString(values[i]));}System.out.print(']');
    }
    private static void flatten4(float[][][][] a) {float[] f=new float[NX*NY*NZ*a[0][0][0].length];int i=0;for(float[][][] p:a)for(float[][] q:p)for(float[] r:q)for(float v:r)f[i++]=v;printArray(f);}
    private static void flatten2(float[][] a) {float[] f=new float[a.length*a[0].length];int i=0;for(float[] r:a)for(float v:r)f[i++]=v;printArray(f);}
    public static void main(String[] args) {
        boolean useField=args.length==0||Boolean.parseBoolean(args[0]);
        float ridge=args.length<2?0.2f:Float.parseFloat(args[1]);
        float[][][] memberships=new float[NX][NY][NZ];
        float[][][][] initial=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)initial[x][y][z][k]=membership(x,y,z,k);
        float[][] variances=imageVariance(useField);float[][] deviations=classStddev(useField);
        float[][] similarity=similarities(variances[0],deviations);float[][][][] psf=bias(useField,ridge);
        float[][][][] data=dataEnergy(useField,variances[0],psf);
        // Run a separate update here so the harness retains the complete output array.
        float[] flat=new float[NX*NY*NZ*CLASSES];int t=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)flat[t++]=initial[x][y][z][k];
        float objective=0.0f;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0.0f;float[] num=new float[CLASSES],energy=new float[CLASSES];boolean[] suppressed=new boolean[CLASSES];
            for(int k=0;k<CLASSES;k++){
                num[k]=data[x][y][z][k];energy[k]=num[k];float ngb=0.0f;int[][] ds={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                for(int[] d:ds)for(int m=0;m<CLASSES;m++)if(m!=k){float v=flat[(((x+d[0])*NY+y+d[1])*NZ+z+d[2])*CLASSES+m];ngb+=v*v;}
                num[k]+=SMOOTH_PRIOR[k]*SMOOTH*ngb;energy[k]+=0.5f*SMOOTH_PRIOR[k]*SMOOTH*ngb;
                float inv=0.0f;for(int m=0;m<CLASSES;m++)if(m!=k)inv+=ATLAS_PRIOR[m]*similarity[k][m]*shape(x,y,z,m)*shape(x,y,z,m);
                num[k]+=PRIOR_COEFFICIENT*inv;energy[k]+=0.5f*PRIOR_COEFFICIENT*inv;
            }
            for(int k=0;k<CLASSES;k++){
                boolean direct=closest(x,y,z,0)==k||closest(x,y,z,1)==k,neighbor=closest(x,y,z,2)==k||closest(x,y,z,3)==k;
                float sh=shape(x,y,z,k),factor;boolean priorOnly=!direct&&!neighbor;
                if(direct)factor=1.0f;else if(neighbor)factor=Math.max(0.25f,Math.min(2.0f*PRIOR_FACTOR*sh,1.0f));else if(sh==0.0f)factor=0.0f;else factor=Math.min(1.0f,2.0f*PRIOR_FACTOR*sh);
                float value;if(priorOnly&&sh==0.0f){value=0.0f;suppressed[k]=true;}else value=factor*(num[k]>ZERO?1.0f/num[k]:INF);
                flat[((x*NY+y)*NZ+z)*CLASSES+k]=value;den+=value;
            }
            if(den>0.0f)for(int k=0;k<CLASSES;k++){int idx=((x*NY+y)*NZ+z)*CLASSES+k;flat[idx]=flat[idx]/den;if(!suppressed[k])objective+=energy[k]*flat[idx]*flat[idx];}
        }
        System.out.print("{\"source\":\"MultipleToadSegmentation.computeExactMemberships optimized tissue one-pass formula composition\",");
        System.out.print("\"segmentation_source_sha256\":\""+SEGMENTATION_SHA+"\",\"use_field\":"+useField+",\"shape\":[7,6,5,3],\"psf_ridge\":"+Float.toString(ridge)+",\"prior_coefficient\":0.15,");
        System.out.print("\"image_variance_weights\":");printArray(variances[0]);System.out.print(",\"image_stddev\":");printArray(variances[1]);
        System.out.print(",\"class_standard_deviations\":");flatten2(deviations);System.out.print(",\"similarity\":[");for(int k=0;k<CLASSES;k++){if(k>0)System.out.print(',');printArray(similarity[k]);}System.out.print(']');
        System.out.print(",\"bias\":");flatten4(psf);System.out.print(",\"data_energy\":");flatten4(data);System.out.print(",\"updated\":");printArray(flat);System.out.println(",\"objective\":"+Float.toString(objective)+"}");
    }
}
