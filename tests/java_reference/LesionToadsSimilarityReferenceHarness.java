/** Selected Java computeMemberships centroid-similarity formula excerpt. */
public final class LesionToadsSimilarityReferenceHarness {
    private static final int CHANNELS=2,CLASSES=4;
    private final String name;
    private final int t1;
    private final float priorScale=0.15f;
    private final float[][] centroid=new float[CHANNELS][CLASSES];
    private final float[][] stddev=new float[CHANNELS][CLASSES];
    private final float[] lesionCentroid={1.15f,2.05f};
    private final float[] lesionStDev={0.4f,0.6f};
    private final float[] channelWeight={0.3f,0.7f};
    private final String[] names={"Sulcal-CSF","Cerebrum-GM","Cerebrum-WM","Cerebellum"};
    private final String[] objectTypes;
    private final float[][] similarity=new float[CLASSES][CLASSES];
    private final float[] lesionSimilarity=new float[CLASSES];
    private final float[] blackHoleSimilarity=new float[CLASSES];

    private LesionToadsSimilarityReferenceHarness(String name) {
        this.name=name;
        t1=name.equals("no_t1")?-1:0;
        objectTypes=name.equals("no_lesion")
            ?new String[]{"obj","obj","obj","obj"}
            :new String[]{"obj","blackhole","lesion","optimized_lesion"};
        for(int c=0;c<CHANNELS;c++)for(int k=0;k<CLASSES;k++) {
            centroid[c][k]=c*1.25f+k*0.7f+((c+k)%2)*0.15f;
            stddev[c][k]=0.55f+0.13f*((c*2+k)%4);
        }
        if(name.equals("zero_std")) {
            stddev[0][0]=0.0f;stddev[0][1]=0.0f;
            centroid[0][1]=centroid[0][0];
            lesionCentroid[0]=centroid[0][0];
            lesionStDev[0]=0.0f;
        }
        int lesionClass=-1;
        for(int k=0;k<CLASSES;k++) {
            lesionSimilarity[k]=1.0f;
            blackHoleSimilarity[k]=1.0f;
            if(objectTypes[k].equals("lesion")||objectTypes[k].equals("blackhole")
                    ||objectTypes[k].equals("optimized_lesion"))lesionClass=k;
        }
        if(lesionClass!=-1) {
            lesionSimilarity[lesionClass]=0.0f;
            blackHoleSimilarity[lesionClass]=0.0f;
        }
    }
    private float similarity(float distance) {
        return 1.0f/(1.0f+(1.0f-priorScale)/priorScale*distance);
    }
    private void calculate() {
        boolean csfFlag=t1!=-1;
        for(int k=0;k<CLASSES;k++)for(int m=0;m<CLASSES;m++) {
            if(m==k){similarity[k][m]=0.0f;continue;}
            float distance=0.0f;
            if(csfFlag&&(names[m].equals("Sulcal-CSF")||names[k].equals("Sulcal-CSF"))) {
                distance+=(centroid[t1][k]-centroid[t1][m])*(centroid[t1][k]-centroid[t1][m])
                    /(stddev[t1][k]+stddev[t1][m])/(stddev[t1][k]+stddev[t1][m]);
            } else {
                for(int c=0;c<CHANNELS;c++)distance+=channelWeight[c]
                    *(centroid[c][k]-centroid[c][m])*(centroid[c][k]-centroid[c][m])
                    /(stddev[c][k]+stddev[c][m])/(stddev[c][k]+stddev[c][m]);
            }
            similarity[k][m]=similarity(distance);
        }
        for(int k=0;k<CLASSES;k++) {
            float distance=0.0f;
            for(int c=0;c<CHANNELS;c++)distance+=channelWeight[c]
                *(lesionCentroid[c]-centroid[c][k])*(lesionCentroid[c]-centroid[c][k])
                /(stddev[c][k]+lesionStDev[c])/(stddev[c][k]+lesionStDev[c]);
            if(!objectTypes[k].equals("lesion")&&!objectTypes[k].equals("blackhole")
                    &&!objectTypes[k].equals("optimized_lesion"))lesionSimilarity[k]=similarity(distance);
            else lesionSimilarity[k]=0.0f;
        }
        // Source blackHoleSimilarity update remains commented out.
    }
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String vector(float[] values){StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(bits(values[i]));}return b.append(']').toString();}
    private static String matrix(float[][] values){float[] out=new float[CLASSES*CLASSES];int q=0;for(float[] row:values)for(float v:row)out[q++]=v;return vector(out);}
    private void emit(boolean first) {
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+name+"\",\"t1_channel\":"+t1+",\"types\":[\""+String.join("\",\"",objectTypes)+"\"],\"centroid\":[");
        for(int c=0;c<CHANNELS;c++){if(c>0)System.out.print(',');System.out.print(vector(centroid[c]));}
        System.out.print("],\"stddev\":[");for(int c=0;c<CHANNELS;c++){if(c>0)System.out.print(',');System.out.print(vector(stddev[c]));}
        System.out.print("],\"lesion_centroid\":"+vector(lesionCentroid)+",\"lesion_stddev\":"+vector(lesionStDev)+",\"channel_weight\":"+vector(channelWeight)+",\"class_similarity\":"+matrix(similarity)+",\"lesion_similarity\":"+vector(lesionSimilarity)+",\"black_hole_similarity\":"+vector(blackHoleSimilarity)+"}");
    }
    public static void main(String[] args) {
        String[] cases={"mixed_t1","no_t1","no_lesion","zero_std"};
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships centroid similarity setup\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        for(int i=0;i<cases.length;i++){LesionToadsSimilarityReferenceHarness h=new LesionToadsSimilarityReferenceHarness(cases[i]);h.calculate();h.emit(i==0);}
        System.out.println("]}");
    }
}
