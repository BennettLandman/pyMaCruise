/** Pinned direct-reference setters and call order from the selected sources. */
public final class LesionToadSegmentationSpecialCentroidSettersReferenceHarness {
    private static final class State {
        float[] lesionCentroid;
        float[] blackHoleCentroid;
        String order = "";
        void setLesionCentroid(float[] cent) { lesionCentroid=cent; order+="lesion,"; }
        void setBlackHoleCentroid(float[] cent) { blackHoleCentroid=cent; order+="black_hole,"; }
    }
    private static String bits(float[] v) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<v.length;i++){if(i>0)s.append(',');s.append(Integer.toUnsignedLong(Float.floatToIntBits(v[i])));}
        return s.append(']').toString();
    }
    public static void main(String[] args) {
        float[] lesion={Float.intBitsToFloat(1036831949),Float.intBitsToFloat(-1072693248),Float.intBitsToFloat(1078530011)};
        float[] blackHole={Float.intBitsToFloat(-2147483648),Float.intBitsToFloat(1051372203),Float.intBitsToFloat(507307272)};
        State state=new State();
        state.setLesionCentroid(lesion);
        state.setBlackHoleCentroid(blackHole);
        System.out.println("{\"source_methods\":[\"LesionToadSegmentation.setLesionCentroid\",\"LesionToadSegmentation.setBlackHoleCentroid\"],\"call_order\":[\"lesion\",\"black_hole\"],\"lesion_input_bits\":"+bits(lesion)+",\"black_hole_input_bits\":"+bits(blackHole)+",\"lesion_same_reference\":"+(state.lesionCentroid==lesion)+",\"black_hole_same_reference\":"+(state.blackHoleCentroid==blackHole)+",\"lesion_bits\":"+bits(state.lesionCentroid)+",\"black_hole_bits\":"+bits(state.blackHoleCentroid)+"}");
    }
}
