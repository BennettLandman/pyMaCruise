/** Selected unfielded OBJECT/LESION/BLACKHOLE energy formula excerpt. */
public final class LesionToadsDataTermsReferenceHarness {
    private static final int CHANNELS=3;
    private static final float[] WEIGHT={0.2f,0.3f,0.5f};
    private static final String[] TYPES={"obj","lesion","blackhole"};
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String vector(float[] values){StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(bits(values[i]));}return b.append(']').toString();}
    private static void run(int typeIndex,int sample,boolean first) {
        float[] image=new float[CHANNELS],classCentroid=new float[CHANNELS];
        float[] lesionCentroid={1.13f,1.91f,2.47f};
        float[] blackHoleCentroid={0.31f,2.63f,1.57f};
        float classEnergy=0.0f,lesionEnergy=0.0f,blackHoleEnergy=0.0f;
        for(int c=0;c<CHANNELS;c++) {
            image[c]=0.85f+sample*0.19f+c*0.73f;
            classCentroid[c]=0.31f*(c+1)+typeIndex*0.17f;
            float difference=image[c]-classCentroid[c];
            classEnergy+=WEIGHT[c]*difference*difference;
            if(typeIndex==1||typeIndex==2) {
                float lesionDifference=image[c]-lesionCentroid[c];
                lesionEnergy+=WEIGHT[c]*lesionDifference*lesionDifference;
            }
            if(typeIndex==2) {
                float blackDifference=image[c]-blackHoleCentroid[c];
                blackHoleEnergy+=WEIGHT[c]*blackDifference*blackDifference;
            }
        }
        if(!first)System.out.print(',');
        System.out.print("{\"type\":\""+TYPES[typeIndex]+"\",\"sample\":"+sample+",\"image\":"+vector(image)+",\"class_centroid\":"+vector(classCentroid)+",\"lesion_centroid\":"+vector(lesionCentroid)+",\"black_hole_centroid\":"+vector(blackHoleCentroid)+",\"channel_weight\":"+vector(WEIGHT)+",\"energies\":"+vector(new float[]{classEnergy,lesionEnergy,blackHoleEnergy})+"}");
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships unfielded OBJECT/LESION/BLACKHOLE data terms\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        boolean first=true;
        for(int type=0;type<TYPES.length;type++)for(int sample=0;sample<3;sample++) {run(type,sample,first);first=false;}
        System.out.println("]}");
    }
}
