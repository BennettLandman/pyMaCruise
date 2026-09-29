import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected LESION/useLesionWeight FLAIR PowerTable data-term excerpt. */
public final class LesionToadsFlairWeightReferenceHarness {
    private static final int CHANNELS=3,FLAIR=1;
    private static final float EXPONENT=0.73f;
    private static final float[] WEIGHT={0.2f,0.3f,0.5f};
    private static final float[] LESION_CENTROID={1.13f,1.91f,2.47f};
    private static final float[] CLASS_CENTROID={0.31f,0.62f,0.93f};
    private static final PowerTable TABLE=new PowerTable(0.0f,4.0f,0.000001f,EXPONENT);
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String vector(float[] values){StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(bits(values[i]));}return b.append(']').toString();}
    private static void run(int sample,boolean useField,boolean first) {
        float[] image=new float[CHANNELS],field=new float[CHANNELS];
        float num=0.0f,numL=0.0f;
        for(int c=0;c<CHANNELS;c++) {
            image[c]=0.75f+sample*1.8f+c*0.61f;
            field[c]=0.7f+c*0.21f+sample*0.13f;
            float value=useField?field[c]*image[c]:image[c];
            float difference=value-CLASS_CENTROID[c];
            num+=WEIGHT[c]*difference*difference;
        }
        float flairValue=useField?field[FLAIR]*image[FLAIR]:image[FLAIR];
        float lesionDifference=flairValue-LESION_CENTROID[FLAIR];
        float squaredResidual=lesionDifference*lesionDifference;
        numL+=TABLE.lookup(squaredResidual,EXPONENT);
        if(!first)System.out.print(',');
        System.out.print("{\"sample\":"+sample+",\"use_field\":"+useField+",\"flair_channel\":"+FLAIR+",\"image\":"+vector(image)+",\"field\":"+vector(field)+",\"class_centroid\":"+vector(CLASS_CENTROID)+",\"lesion_centroid\":"+vector(LESION_CENTROID)+",\"channel_weight\":"+vector(WEIGHT)+",\"exponent_bits\":"+bits(EXPONENT)+",\"squared_residual_bits\":"+bits(squaredResidual)+",\"energies\":"+vector(new float[]{num,numL,0.0f})+"}");
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships LESION FLAIR useLesionWeight branch + PowerTable\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"powertable_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"cases\":[");
        boolean first=true;
        for(int sample=0;sample<3;sample++)for(boolean useField:new boolean[]{false,true}){run(sample,useField,first);first=false;}
        System.out.println("]}");
    }
}
