/** Selected computeExactMemberships OUTLIER smoothing/prior/inversion excerpt. */
public final class ToadsOutlierMembershipUpdateReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float[][] CASES={
        {0.7f,0.9f,0.2f,1.0f,0.4f,0.6f,2.5f,0.3f,1.7f},
        {1.25f,0.35f,0.65f,0.25f,0.125f,0.8f,1.75f,0.2f,2.25f},
        {0.5f,1.1f,0.875f,0.8125f,0.7f,0.15f,3.125f,0.33f,0.625f},
        {0.0f,0.0f,0.5f,1.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
        {-0.1f,0.4f,0.25f,0.5f,0.2f,0.3f,0.75f,0.1f,0.5f},
        {0.3f,0.8f,0.4f,0.625f,0.0f,0.0f,0.0f,0.0f,0.7f}
    };
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static float[] run(float[] a){
        float num=a[0],numO=a[1],outlier=a[2],factor=a[3];
        float energy=(1.0f-outlier)*num+outlier*numO;
        float smoothPrior=a[4],smoothing=a[5],ngb=a[6],priorCoefficient=a[7],invprior=a[8];
        num+=smoothPrior*smoothing*ngb;numO+=smoothPrior*smoothing*ngb;
        energy+=0.5f*smoothPrior*smoothing*ngb;
        num+=priorCoefficient*invprior;numO+=priorCoefficient*invprior;
        energy+=0.5f*priorCoefficient*invprior;
        if(num>1.0e-10f)num=1.0f/num;else num=1.0e10f;
        if(numO>1.0e-10f)numO=1.0f/numO;else numO=1.0e10f;
        float mem=factor*((1.0f-outlier)*num+outlier*numO);
        float updatedOutlier=factor*outlier*numO;
        return new float[]{mem,updatedOutlier,energy};
    }
    public static void main(String[] args){
        StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        for(int i=0;i<CASES.length;i++){
            if(i>0)s.append(',');float[] a=CASES[i],r=run(a);
            s.append("{\"input_bits\":[");for(int j=0;j<a.length;j++){if(j>0)s.append(',');s.append(bits(a[j]));}
            s.append("],\"class_bits\":").append(bits(r[0])).append(",\"outlier_bits\":").append(bits(r[1]))
                .append(",\"energy_bits\":").append(bits(r[2])).append('}');
        }
        System.out.println(s.append("]}").toString());
    }
}
