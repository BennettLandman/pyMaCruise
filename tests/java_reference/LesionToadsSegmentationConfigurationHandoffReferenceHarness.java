/** Pinned configuration setters and call order from selected Java sources. */
public final class LesionToadsSegmentationConfigurationHandoffReferenceHarness {
    private static final class State {
        boolean useLesionWeight=false, gmFlag=false;
        float[] intensityMax={1.0f,1.0f}, intensityScale={1.0f,1.0f};
        void setLesionWeight(){useLesionWeight=true;}
        void reSetLesionWeight(){useLesionWeight=false;}
        void setIntensityMax(float[] values){intensityMax=values;}
        void setIntensityScale(float[] values){intensityScale=values;}
        void setGMflag(){gmFlag=true;}
    }
    private static final float[] MAX={Float.intBitsToFloat(1065353216),Float.intBitsToFloat(1075838976)};
    private static final float[] SCALE={Float.intBitsToFloat(1056964608),Float.intBitsToFloat(1082130432)};
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static String vector(float[] a){return "["+Integer.toUnsignedLong(bits(a[0]))+","+Integer.toUnsignedLong(bits(a[1]))+"]";}
    private static String scenario(boolean lesionWeight){
        State state=new State();
        if(lesionWeight)state.setLesionWeight();else state.reSetLesionWeight();
        state.setIntensityMax(MAX);
        state.setIntensityScale(SCALE);
        state.setGMflag();
        return "{\"lesion_weight_input\":"+lesionWeight+",\"use_lesion_weight\":"+state.useLesionWeight+",\"gm_flag\":"+state.gmFlag+",\"intensity_max_same_reference\":"+(state.intensityMax==MAX)+",\"intensity_scale_same_reference\":"+(state.intensityScale==SCALE)+",\"intensity_max_bits\":"+vector(state.intensityMax)+",\"intensity_scale_bits\":"+vector(state.intensityScale)+"}";
    }
    public static void main(String[] args){
        System.out.println("{\"source_methods\":[\"setLesionWeight/reSetLesionWeight\",\"setIntensityMax\",\"setIntensityScale\",\"setGMflag\"],\"call_order\":[\"lesion_weight_branch\",\"setIntensityMax\",\"setIntensityScale\",\"setGMflag\"],\"max_input_bits\":"+vector(MAX)+",\"scale_input_bits\":"+vector(SCALE)+",\"false_case\":"+scenario(false)+",\"true_case\":"+scenario(true)+"}");
    }
}
