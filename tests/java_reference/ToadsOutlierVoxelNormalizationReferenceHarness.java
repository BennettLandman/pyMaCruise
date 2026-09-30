/** Selected computeExactMemberships final denominator/energy excerpt. */
public final class ToadsOutlierVoxelNormalizationReferenceHarness {
    private static final String SOURCE_SHA =
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String[][] TYPES={{"OBJECT","OUTLIER","MASK"},{"OUTLIER","OBJECT","OUTLIER"},{"OBJECT","MASK"}};
    private static final float[][] MEMBERS={{0.13f,0.37f,0.22f},{0.41f,0.17f,0.29f},{0.0f,0.0f}};
    private static final float[][] ENERGIES={{0.7f,1.1f,0.35f},{0.25f,0.8f,1.25f},{1.0f,2.0f}};
    private static final float[] OUTLIERS={0.19f,0.31f,0.5f};
    private static int bits(float f){return Float.floatToIntBits(f);}
    private static void ints(StringBuilder s,float[] a){
        s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}s.append(']');
    }
    private static void strings(StringBuilder s,String[] a){
        s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append('"').append(a[i]).append('"');}s.append(']');
    }
    private static float[] run(float[] m,float[] e,String[] types,float outlier){
        float den=0.0f;for(int k=0;k<m.length;k++)den+=m[k];
        float[] result=new float[m.length+3];result[0]=den;float dist=0.0f;
        if(den>0.0f){
            for(int k=0;k<m.length;k++){
                m[k]=m[k]/den;
                if(types[k].equals("OUTLIER")||types[k].equals("OUTMASK"))outlier=outlier/den;
                float t=e[k]*m[k];t=t*m[k];dist+=t;
            }
        }
        for(int k=0;k<m.length;k++)result[k+1]=m[k];result[m.length+1]=outlier;result[m.length+2]=dist;return result;
    }
    public static void main(String[] args){
        StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"cases\":[");
        for(int i=0;i<TYPES.length;i++){
            if(i>0)s.append(',');float[] mem=MEMBERS[i].clone(),en=ENERGIES[i].clone();float out=OUTLIERS[i];
            float[] outputs=run(mem,en,TYPES[i],out);
            s.append("{\"class_types\":");strings(s,TYPES[i]);s.append(",\"membership_bits\":");ints(s,MEMBERS[i]);
            s.append(",\"energy_bits\":");ints(s,ENERGIES[i]);s.append(",\"outlier_bits\":").append(bits(out))
                .append(",\"denominator_bits\":").append(bits(outputs[0])).append(",\"normalized_bits\":");
            float[] normalized=new float[TYPES[i].length];System.arraycopy(outputs,1,normalized,0,normalized.length);ints(s,normalized);
            s.append(",\"normalized_outlier_bits\":").append(bits(outputs[TYPES[i].length+1]))
                .append(",\"objective_bits\":").append(bits(outputs[TYPES[i].length+2])).append('}');
        }
        System.out.println(s.append("]}").toString());
    }
}
