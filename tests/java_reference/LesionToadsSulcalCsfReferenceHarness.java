import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected Sulcal-CSF/MPRAGE ordinary-object data-energy excerpt. */
public final class LesionToadsSulcalCsfReferenceHarness {
    private static final int CHANNELS=3,T1=0;
    private static final float EXPONENT=0.61f;
    private static final float[] WEIGHT={0.2f,0.3f,0.5f};
    private static final PowerTable TABLE=new PowerTable(0.0f,4.0f,0.000001f,EXPONENT);
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String vector(float[] values){StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(bits(values[i]));}return b.append(']').toString();}
    private static void run(String name,int sample,boolean useField,boolean hasMprage,int t1,boolean first) {
        float[] image=new float[CHANNELS],field=new float[CHANNELS],centroid=new float[CHANNELS];
        float energy=0.0f;
        for(int c=0;c<CHANNELS;c++) {
            image[c]=0.75f+sample*1.8f+c*0.61f;
            field[c]=0.7f+c*0.21f+sample*0.13f;
            centroid[c]=0.31f*(c+1);
        }
        if(t1!=-1&&name.equals("Sulcal-CSF")) {
            float value=useField?field[t1]*image[t1]:image[t1];
            float difference=value-centroid[t1];
            float squared=difference*difference;
            energy=hasMprage?(float)TABLE.lookup(squared,EXPONENT):squared;
        } else {
            for(int c=0;c<CHANNELS;c++) {
                float value=useField?field[c]*image[c]:image[c];
                float difference=value-centroid[c];
                energy+=WEIGHT[c]*difference*difference;
            }
        }
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+name+"\",\"sample\":"+sample+",\"use_field\":"+useField+",\"has_mprage\":"+hasMprage+",\"t1_channel\":"+t1+",\"image\":"+vector(image)+",\"field\":"+vector(field)+",\"class_centroid\":"+vector(centroid)+",\"channel_weight\":"+vector(WEIGHT)+",\"exponent_bits\":"+bits(EXPONENT)+",\"energy_bits\":"+bits(energy)+"}");
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships Sulcal-CSF/MPRAGE data branch + PowerTable\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"powertable_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"cases\":[");
        Object[][] variants={{"Sulcal-CSF",0,false,true,0},{"Sulcal-CSF",1,false,true,0},{"Sulcal-CSF",2,false,true,0},{"Sulcal-CSF",0,true,true,0},{"Sulcal-CSF",2,true,true,0},{"Sulcal-CSF",1,false,false,0},{"Sulcal-CSF",1,true,false,0},{"Sulcal-CSF",1,false,false,-1},{"Cerebrum-GM",1,false,true,0}};
        for(int i=0;i<variants.length;i++){Object[] v=variants[i];run((String)v[0],(int)v[1],(boolean)v[2],(boolean)v[3],(int)v[4],i==0);}
        System.out.println("]}");
    }
}
