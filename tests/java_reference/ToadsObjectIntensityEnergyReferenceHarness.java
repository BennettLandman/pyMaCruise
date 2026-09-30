/** Selected computeExactMemberships OBJECT residual-energy excerpt. */
public final class ToadsObjectIntensityEnergyReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX = 4, NY = 3, NZ = 2, CHANNELS = 3, CLASSES = 4;
    private static int code(int index, int modulus, int offset) { return (index * 7 + offset) % modulus; }
    private static float[][][][] images() {
        float[][][][] a = new float[NX][NY][NZ][CHANNELS]; int n = 0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++,n++)
            for (int c=0;c<CHANNELS;c++) a[x][y][z][c] = (code(n+c,17,5)-8)/9.0f;
        return a;
    }
    private static float[][] centroids() {
        float[][] a = new float[CHANNELS][CLASSES];
        for (int c=0;c<CHANNELS;c++) for (int k=0;k<CLASSES;k++) a[c][k] = ((c*11+k*3)%13-6)/10.0f;
        return a;
    }
    private static float[] modalityWeights() { return new float[] {0.7f, 1.1f, 0.35f}; }
    private static float[][] classWeights() {
        float[][] a = new float[CHANNELS][CLASSES];
        for (int c=0;c<CHANNELS;c++) for (int k=0;k<CLASSES;k++) a[c][k] = 0.5f+(c+k)/10.0f;
        return a;
    }
    private static float[] imageVariances() { return new float[] {0.8f, 1.25f, 0.625f}; }
    private static float[][][][] field() {
        float[][][][] a = new float[NX][NY][NZ][CHANNELS]; int n=0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++,n++)
            for (int c=0;c<CHANNELS;c++) a[x][y][z][c] = 0.8f+((n+c)%5)*0.1f;
        return a;
    }
    private static float[][][][] bias() {
        float[][][][] a = new float[NX][NY][NZ][CHANNELS]; int n=0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++,n++)
            for (int c=0;c<CHANNELS;c++) a[x][y][z][c] = ((n*3+c)%9)/23.0f;
        return a;
    }
    private static float[][][][] energy(float[][][][] image, float[][] means, float[] modality,
            float[][] classVar, float[] imageVar, float[][][][] field, float[][][][] bias) {
        float[][][][] result = new float[NX][NY][NZ][CLASSES];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            for (int k=0;k<CLASSES;k++) for (int c=0;c<CHANNELS;c++) {
                float corrected = field == null ? image[x][y][z][c] : field[x][y][z][c]*image[x][y][z][c];
                float residual = corrected-(bias == null ? 0.0f : bias[x][y][z][c]);
                residual = residual-means[c][k];
                float weight = modality[c]*classVar[c][k]; weight = weight*imageVar[c];
                float contribution = weight*residual; contribution = contribution*residual;
                result[x][y][z][k] = result[x][y][z][k]+contribution;
            }
        return result;
    }
    private static void appendBits(StringBuilder out, float[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if(i>0)out.append(',');out.append(Float.floatToIntBits(values[i])); }
        out.append(']');
    }
    private static void append3D(StringBuilder out, float[][][][] a, int finalAxis) {
        float[] values = new float[NX*NY*NZ*finalAxis]; int p=0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) for (int c=0;c<finalAxis;c++) values[p++]=a[x][y][z][c];
        appendBits(out, values);
    }
    private static void append2D(StringBuilder out, float[][] a) {
        float[] values=new float[a.length*a[0].length];int p=0;for(float[] row:a)for(float v:row)values[p++]=v;appendBits(out,values);
    }
    private static void append1D(StringBuilder out, float[] a) { appendBits(out,a); }
    private static void emit(StringBuilder out, String name, boolean useField, boolean useBias) {
        float[][][][] image=images(), f=useField?field():null, b=useBias?bias():null;
        float[][] means=centroids(), classVar=classWeights();float[] modality=modalityWeights(),imageVar=imageVariances();
        float[][][][] result=energy(image,means,modality,classVar,imageVar,f,b);
        if(out.charAt(out.length()-1)!='[')out.append(',');
        out.append("{\"name\":\"").append(name).append("\",\"shape\":[").append(NX).append(',').append(NY).append(',').append(NZ).append(',').append(CHANNELS)
            .append("],\"classes\":").append(CLASSES).append(",\"images\":");append3D(out,image,CHANNELS);
        out.append(",\"centroids\":");append2D(out,means);out.append(",\"modality_weights\":");append1D(out,modality);
        out.append(",\"class_variance_weights\":");append2D(out,classVar);out.append(",\"image_variances\":");append1D(out,imageVar);
        out.append(",\"field\":");if(f==null)out.append("null");else append3D(out,f,CHANNELS);
        out.append(",\"psf_bias\":");if(b==null)out.append("null");else append3D(out,b,CHANNELS);
        out.append(",\"energy_bits\":");append3D(out,result,CLASSES);out.append('}');
    }
    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        emit(out,"raw_no_bias",false,false);emit(out,"field_no_bias",true,false);
        emit(out,"raw_psf_bias",false,true);emit(out,"field_psf_bias",true,true);
        System.out.println(out.append("]}"));
    }
}
