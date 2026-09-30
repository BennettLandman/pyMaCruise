/** Selected Java source-formula reference for a composed OBJECT/MASK/OUTLIER pass. */
public final class ToadsObjectMaskOutlierMembershipPassReferenceHarness {
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=4,NY=4,NZ=4,CLASSES=3,OUT=2;
    private static final float SMOOTH=0.14f,FACTOR=0.58f,COEFFICIENT=0.22f,OUTLIER_VALUE=1.3f,ZERO=1e-10f,INF=1e10f;
    private static final float MODALITY=0.9f,IMAGE_VARIANCE=1.2f;
    private static final float[] VARIANCE={0.8f,1.1f,0.7f},MEAN={-0.25f,0.15f,0.4f};
    private static final float[] SMOOTH_PRIOR={0.6f,0.9f,1.2f},ATLAS_PRIOR={0.2f,0.3f,0.5f};
    private static final float[][] SIM={{0,.4f,.2f},{.3f,0,.5f},{.25f,.45f,0}};
    private static float image(int x,int y,int z){return (((x+2*y+3*z)%13)-6)/7.0f;}
    private static float field(int x,int y,int z){return 0.8f+((2*x+y+z)%5)*0.05f;}
    private static float bias(int x,int y,int z){return ((x+3*y+2*z)%9)/23.0f;}
    private static float shape(int x,int y,int z,int k){return 0.3f+((2*x+3*y+z+k)%11)*0.04f;}
    private static float initial(int x,int y,int z,int k){return 0.1f+((x+2*y+3*z+k)%9)*0.02f;}
    private static float outlier(int x,int y,int z){return 0.15f+((x+2*y+z)%10)*0.025f;}
    private static float weight(int k){return (MODALITY*VARIANCE[k])*IMAGE_VARIANCE;}
    private static float[][][][] update(float[][][] outMap,float[] objectiveOut){
        float[][][][] memberships=new float[NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)memberships[x][y][z][k]=initial(x,y,z,k);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)outMap[x][y][z]=outlier(x,y,z);
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        float total=0;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0;float[] energy=new float[CLASSES];
            for(int k=0;k<CLASSES;k++){
                float value=image(x,y,z),corrected=field(x,y,z)*value-bias(x,y,z),residual=corrected-MEAN[k];
                float num;
                if(k==1)num=weight(k)*(value==0.0f?0.0f:1.0f);
                else num=weight(k)*residual*residual;
                float numO=0;
                if(k==OUT){numO=weight(k)*OUTLIER_VALUE;energy[k]=(1.0f-outMap[x][y][z])*num+outMap[x][y][z]*numO;}
                else energy[k]=num;
                float objective=energy[k];
                float neighbor=0;
                for(int[] d:offsets)for(int m=0;m<CLASSES;m++)if(m!=k){float v=memberships[x+d[0]][y+d[1]][z+d[2]][m];neighbor=neighbor+v*v;}
                float add=SMOOTH_PRIOR[k]*SMOOTH*neighbor;
                num=num+add;if(k==OUT)numO=numO+add;
                objective=objective+(0.5f*SMOOTH_PRIOR[k])*SMOOTH*neighbor;
                float invprior=0;
                for(int m=0;m<CLASSES;m++)if(m!=k){float term=ATLAS_PRIOR[m]*SIM[k][m];term=term*shape(x,y,z,m);term=term*shape(x,y,z,m);invprior=invprior+term;}
                num=num+COEFFICIENT*invprior;if(k==OUT)numO=numO+COEFFICIENT*invprior;
                objective=objective+(0.5f*COEFFICIENT)*invprior;
                float shapeValue=shape(x,y,z,k),factor=Math.min(1.0f,2.0f*FACTOR*shapeValue);
                float reciprocal=num>ZERO?1.0f/num:INF,mem;
                if(k==OUT){float reciprocalO=numO>ZERO?1.0f/numO:INF;mem=factor*((1.0f-outMap[x][y][z])*reciprocal+outMap[x][y][z]*reciprocalO);outMap[x][y][z]=(factor*outMap[x][y][z])*reciprocalO;}
                else mem=factor*reciprocal;
                memberships[x][y][z][k]=mem;den=den+mem;energy[k]=objective;
            }
            if(den>0)for(int k=0;k<CLASSES;k++){memberships[x][y][z][k]=memberships[x][y][z][k]/den;if(k==OUT)outMap[x][y][z]=outMap[x][y][z]/den;total=total+energy[k]*memberships[x][y][z][k]*memberships[x][y][z][k];}
        }
        objectiveOut[0]=total;return memberships;
    }
    private static void append(StringBuilder o,int[] a){o.append('[');for(int i=0;i<a.length;i++){if(i>0)o.append(',');o.append(a[i]);}o.append(']');}
    private static int[] bits(float[][][][] a){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(a[x][y][z][k]);return b;}
    private static int[] bits(float[][][] a){int[] b=new int[NX*NY*NZ];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)b[p++]=Float.floatToIntBits(a[x][y][z]);return b;}
    private static int[] volume(int mode){int[] b=new int[NX*NY*NZ];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){float v=mode==0?image(x,y,z):mode==1?field(x,y,z):mode==2?bias(x,y,z):outlier(x,y,z);b[p++]=Float.floatToIntBits(v);}return b;}
    private static int[] shapes(){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(shape(x,y,z,k));return b;}
    private static int[] initial(){int[] b=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)b[p++]=Float.floatToIntBits(initial(x,y,z,k));return b;}
    public static void main(String[] args){float[][][] out=new float[NX][NY][NZ];float[] objective={0};float[][][][] updated=update(out,objective);StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[4,4,4],\"classes\":3,\"class_types\":[\"OBJECT\",\"MASK\",\"OUTLIER\"]");
        s.append(",\"images_bits\":");append(s,volume(0));s.append(",\"field_bits\":");append(s,volume(1));s.append(",\"bias_bits\":");append(s,volume(2));s.append(",\"initial_outlier_bits\":");append(s,volume(3));s.append(",\"initial_memberships_bits\":");append(s,initial());s.append(",\"shape_bits\":");append(s,shapes());
        s.append(",\"centroids_bits\":");append(s,new int[]{Float.floatToIntBits(MEAN[0]),Float.floatToIntBits(MEAN[1]),Float.floatToIntBits(MEAN[2])});s.append(",\"modality_weight_bits\":").append(Float.floatToIntBits(MODALITY));s.append(",\"variance_bits\":");append(s,new int[]{Float.floatToIntBits(VARIANCE[0]),Float.floatToIntBits(VARIANCE[1]),Float.floatToIntBits(VARIANCE[2])});s.append(",\"image_variance_bits\":").append(Float.floatToIntBits(IMAGE_VARIANCE));
        s.append(",\"smoothing_priors_bits\":");append(s,new int[]{Float.floatToIntBits(SMOOTH_PRIOR[0]),Float.floatToIntBits(SMOOTH_PRIOR[1]),Float.floatToIntBits(SMOOTH_PRIOR[2])});s.append(",\"atlas_priors_bits\":");append(s,new int[]{Float.floatToIntBits(ATLAS_PRIOR[0]),Float.floatToIntBits(ATLAS_PRIOR[1]),Float.floatToIntBits(ATLAS_PRIOR[2])});
        int[] sim=new int[9];int p=0;for(float[] row:SIM)for(float v:row)sim[p++]=Float.floatToIntBits(v);s.append(",\"similarity_bits\":");append(s,sim);
        int[] labels=new int[NX*NY*NZ*4];java.util.Arrays.fill(labels,-1);s.append(",\"closest_labels\":");append(s,labels);
        s.append(",\"outlier_value\":").append(OUTLIER_VALUE).append(",\"smoothing\":").append(SMOOTH).append(",\"prior_factor\":").append(FACTOR).append(",\"prior_coefficient\":").append(COEFFICIENT);
        s.append(",\"updated_bits\":");append(s,bits(updated));s.append(",\"updated_outlier_bits\":");append(s,bits(out));s.append(",\"objective_bits\":").append(Float.floatToIntBits(objective[0])).append('}');System.out.println(s);
    }
}
