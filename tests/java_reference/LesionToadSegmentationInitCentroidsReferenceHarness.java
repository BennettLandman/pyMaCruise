/** Pinned behavior of LesionToadSegmentation.initCentroids(float[][]). */
public final class LesionToadSegmentationInitCentroidsReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }
    private static String bitsJson(float[][] values) {
        StringBuilder out=new StringBuilder("[");
        for(int c=0;c<values.length;c++){
            if(c>0)out.append(','); out.append('[');
            for(int k=0;k<values[c].length;k++){
                if(k>0)out.append(','); out.append(Integer.toUnsignedLong(bits(values[c][k])));
            }
            out.append(']');
        }
        return out.append(']').toString();
    }
    public static void main(String[] args) {
        int nc=2, classes=4;
        float[][] input={{0.1f,-2.25f,3.1415927f,0.0f},{-0.0f,0.33333334f,1.0e-20f,16777216.0f}};
        float[][] centroid=new float[nc][classes];
        float[][] priorCentroid=new float[nc][classes];
        for(int c=0;c<nc;c++)for(int k=0;k<classes;k++){
            centroid[c][k]=input[c][k];
            priorCentroid[c][k]=input[c][k];
        }
        System.out.println("{\"source_method\":\"LesionToadSegmentation.initCentroids\",\"shape\":[2,4],\"input_bits\":"+bitsJson(input)+",\"centroid_bits\":"+bitsJson(centroid)+",\"prior_centroid_bits\":"+bitsJson(priorCentroid)+"}");
    }
}
