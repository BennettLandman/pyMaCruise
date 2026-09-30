/** Selected computeExactMemberships OUTLIER data-energy excerpt. */
public final class ToadsOutlierIntensityEnergyReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=3, NY=2, NZ=2, CHANNELS=3, CLASSES=3;
    private static float[][][][] image() {
        float[][][][] a=new float[NX][NY][NZ][CHANNELS];int n=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++,n++)
            for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=((n*5+c*3)%17-8)/11.0f;
        return a;
    }
    private static float[][] means() {
        float[][] a=new float[CHANNELS][CLASSES];
        for(int c=0;c<CHANNELS;c++)for(int k=0;k<CLASSES;k++)a[c][k]=((c*7+k*4)%13-6)/9.0f;
        return a;
    }
    private static float[] modality(){return new float[]{0.7f,1.1f,0.35f};}
    private static float[][] classWeight(){
        float[][] a=new float[CHANNELS][CLASSES];
        for(int c=0;c<CHANNELS;c++)for(int k=0;k<CLASSES;k++)a[c][k]=0.5f+(c+2*k)/12.0f;
        return a;
    }
    private static float[] imageVariance(){return new float[]{0.8f,1.25f,0.625f};}
    private static float[][][][] field(){
        float[][][][] a=new float[NX][NY][NZ][CHANNELS];int n=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++,n++)
            for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=0.75f+((n+2*c)%6)*0.08f;
        return a;
    }
    private static float[][][][] bias(){
        float[][][][] a=new float[NX][NY][NZ][CHANNELS];int n=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++,n++)
            for(int c=0;c<CHANNELS;c++)a[x][y][z][c]=((3*n+c)%11)/29.0f;
        return a;
    }
    private static float[][][][][] energies(float[][][][] im,float[][] centroid,float[] mod,
            float[][] var,float[] imageVar,float[][][][] f,float[][][][] b,float outlier) {
        float[][][][][] result=new float[2][NX][NY][NZ][CLASSES];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)
            for(int k=0;k<CLASSES;k++)for(int c=0;c<CHANNELS;c++) {
                float corrected=f==null?im[x][y][z][c]:f[x][y][z][c]*im[x][y][z][c];
                float residual=corrected-(b==null?0.0f:b[x][y][z][c]);
                residual=residual-centroid[c][k];
                float weight=mod[c]*var[c][k];weight=weight*imageVar[c];
                float contribution=weight*residual;contribution=contribution*residual;
                result[0][x][y][z][k]=result[0][x][y][z][k]+contribution;
                float outContribution=weight*outlier;
                result[1][x][y][z][k]=result[1][x][y][z][k]+outContribution;
            }
        return result;
    }
    private static void bits(StringBuilder s,float[] a){
        s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(a[i]));}s.append(']');
    }
    private static void bits4(StringBuilder s,float[][][][] a){
        float[] v=new float[NX*NY*NZ*CHANNELS];int p=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)v[p++]=a[x][y][z][c];bits(s,v);
    }
    private static void bits2(StringBuilder s,float[][] a){float[] v=new float[a.length*a[0].length];int p=0;for(float[] r:a)for(float f:r)v[p++]=f;bits(s,v);}
    private static void bits3(StringBuilder s,float[][][][] a){
        float[] v=new float[NX*NY*NZ*CLASSES];int p=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)v[p++]=a[x][y][z][k];bits(s,v);
    }
    private static void appendCase(StringBuilder s,String name,boolean useField,boolean useBias,float outlier){
        float[][][][] im=image(),f=useField?field():null,b=useBias?bias():null;
        float[][] m=means(),v=classWeight();float[] mod=modality(),iv=imageVariance();
        float[][][][][] e=energies(im,m,mod,v,iv,f,b,outlier);
        if(s.charAt(s.length()-1)!='[')s.append(',');
        s.append("{\"name\":\"").append(name).append("\",\"outlier\":").append(outlier)
            .append(",\"images_bits\":");bits4(s,im);s.append(",\"centroids_bits\":");bits2(s,m);
        s.append(",\"modality_bits\":");bits(s,mod);s.append(",\"class_weight_bits\":");bits2(s,v);
        s.append(",\"image_variance_bits\":");bits(s,iv);s.append(",\"field_bits\":");if(f==null)s.append("null");else bits4(s,f);
        s.append(",\"bias_bits\":");if(b==null)s.append("null");else bits4(s,b);
        s.append(",\"residual_bits\":");bits3(s,e[0]);s.append(",\"outlier_bits\":");bits3(s,e[1]);s.append('}');
    }
    public static void main(String[] args){
        StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[")
            .append(NX).append(',').append(NY).append(',').append(NZ).append(',').append(CHANNELS)
            .append("],\"classes\":").append(CLASSES).append(",\"cases\":[");
        appendCase(s,"raw",false,false,0.37f);appendCase(s,"field",true,false,0.37f);
        appendCase(s,"bias",false,true,1.125f);appendCase(s,"field_bias",true,true,1.125f);
        System.out.println(s.append("]}").toString());
    }
}
