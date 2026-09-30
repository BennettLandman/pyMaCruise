/** Selected computeExactMemberships OUTLIER relation mixture excerpt. */
public final class ToadsOutlierMembershipMixReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static float[][] cases(){return new float[][]{
        {0.3125f,0.1875f,0.6875f,1.0f},
        {0.8125f,0.625f,0.4375f,0.25f},
        {0.15625f,0.875f,0.53125f,0.8125f},
        {0.9375f,0.3125f,0.21875f,0.33333334f},
        {0.001f,0.999f,0.001f,1.0f},
        {0.5f,0.5f,0.5f,0.5f}
    };}
    private static int bits(float f){return Float.floatToIntBits(f);}
    public static void main(String[] args){
        float[][] values=cases();StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        for(int i=0;i<values.length;i++){
            if(i>0)s.append(',');float n=values[i][0],m=values[i][1],no=values[i][2],factor=values[i][3];
            float mixed=(1.0f-m)*n+m*no;float classResult=factor*mixed;
            float factoredOutlier=factor*m;float outlierResult=factoredOutlier*no;
            s.append("{\"input_bits\":[").append(bits(n)).append(',').append(bits(m)).append(',').append(bits(no)).append(',').append(bits(factor))
                .append("],\"class_bits\":").append(bits(classResult)).append(",\"outlier_bits\":").append(bits(outlierResult)).append('}');
        }
        System.out.println(s.append("]}").toString());
    }
}
