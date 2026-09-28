import java.util.Random;
import edu.jhmi.rad.medic.utilities.PowerTable;

/** Compare selected TOADS optimized-likelihood powers with source PowerTable. */
public final class ToadsPowerTableReferenceHarness {
    private static final String HASH="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final float[] EXPONENTS={0.80f,1.00f,1.00f};

    private static String floats(float[] values) {
        StringBuilder result=new StringBuilder("[");
        for(int i=0;i<values.length;i++) { if(i>0)result.append(','); result.append(Float.toString(values[i])); }
        return result.append(']').toString();
    }

    private static String doubles(double[] values) {
        StringBuilder result=new StringBuilder("[");
        for(int i=0;i<values.length;i++) { if(i>0)result.append(','); result.append(Double.toString(values[i])); }
        return result.append(']').toString();
    }

    private static float[] membership(float[] p, PowerTable[] tables) {
        float[] inverse=new float[3]; float den=0.0f;
        for(int k=0;k<3;k++) {
            float energy=(float)tables[k].lookup((double)p[k],(double)EXPONENTS[k]);
            inverse[k]=energy>1e-10f ? 1.0f/energy : 1e10f;
            den+=inverse[k];
        }
        for(int k=0;k<3;k++) inverse[k]/=den;
        return inverse;
    }

    public static void main(String[] args) {
        PowerTable[] tables=new PowerTable[3];
        for(int k=0;k<3;k++) tables[k]=new PowerTable(0.0f,4.0f,0.000001f,EXPONENTS[k]);
        float[] edge={0.0f,1e-12f,1e-8f,0.49e-6f,0.5e-6f,0.51e-6f,
                1e-6f,1.49e-6f,1.5e-6f,1.51e-6f,0.0001f,0.01f,
                0.25f,0.5f,1.0f,2.0f,3.9999995f,4.0f};
        Random random=new Random(1804L);
        int randomCount=1024;
        int count=edge.length+randomCount;
        System.out.print("{\"source\":\"VandyCRUISE PowerTable.lookup\",\"source_sha256\":\"");
        System.out.print(HASH+"\",\"input_range\":[0,4],\"step\":"+Double.toString((double)0.000001f)+",\"seed\":1804,\"exponents\":[0.8,1,1],\"cases\":[");
        for(int i=0;i<count;i++) {
            float p0=i<edge.length?edge[i]:4.0f*random.nextFloat();
            float p1=(p0*0.61803395f)%4.0f;
            float p2=4.0f-p0;
            float[] p={p0,p1,p2}; double[] powered=new double[3];
            for(int k=0;k<3;k++) powered[k]=tables[k].lookup((double)p[k],(double)EXPONENTS[k]);
            float[] membership=membership(p,tables);
            if(i>0)System.out.print(',');
            System.out.print("{\"p\":"+floats(p)+",\"power\":"+doubles(powered)+",\"membership\":"+floats(membership)+"}");
        }
        System.out.println("]}");
    }
}
