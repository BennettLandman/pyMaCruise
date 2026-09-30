/** Selected computeExactMemberships MASK / maskDistance formula excerpt. */
public final class ToadsMaskIntensityEnergyReferenceHarness {
    private static final String SOURCE_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=3,NY=2,NZ=4,CHANNELS=3,CLASSES=4;
    private static final float[] INTENSITIES={0.0f,-0.0f,1.0f,-1.5f,2.25f,0.01f,0.0f,-2.0f};
    private static final float[] MODALITY={0.75f,1.25f,0.4f};
    private static final float[] VARIANCE={0.6f,1.1f,0.9f,0.3f,1.2f,0.7f,1.05f,0.85f,0.95f,0.5f,0.8f,1.4f};
    private static final float[] IMAGE_VARIANCE={0.8f,1.25f,0.625f};
    private static float image(int x,int y,int z,int c){int n=(x*NY+y)*NZ+z;return INTENSITIES[(n*3+c*5)%INTENSITIES.length];}
    private static float distance(float image){return image==0.0f?0.0f:1.0f;}
    private static int[] flattenImageBits(){int[] a=new int[NX*NY*NZ*CHANNELS];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)a[p++]=Float.floatToIntBits(image(x,y,z,c));return a;}
    private static int[] flattenEnergyBits(){int[] a=new int[NX*NY*NZ*CLASSES];int p=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++){
        float num=0.0f;for(int c=0;c<CHANNELS;c++){float term=MODALITY[c]*VARIANCE[c*CLASSES+k];term=term*IMAGE_VARIANCE[c];term=term*distance(image(x,y,z,c));num=num+term;}a[p++]=Float.floatToIntBits(num);
    }return a;}
    private static int[] floatBits(float[] values){int[] a=new int[values.length];for(int i=0;i<a.length;i++)a[i]=Float.floatToIntBits(values[i]);return a;}
    private static void append(StringBuilder out,int[] a){out.append('[');for(int i=0;i<a.length;i++){if(i>0)out.append(',');out.append(a[i]);}out.append(']');}
    public static void main(String[] args){StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"shape\":[3,2,4,3],\"classes\":4,\"images_bits\":");append(out,flattenImageBits());out.append(",\"modality_weights_bits\":");append(out,floatBits(MODALITY));out.append(",\"class_variance_bits\":");append(out,floatBits(VARIANCE));out.append(",\"image_variances_bits\":");append(out,floatBits(IMAGE_VARIANCE));out.append(",\"energy_bits\":");append(out,flattenEnergyBits());System.out.println(out.append('}'));}
}
