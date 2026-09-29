/**
 * Selected-constructor formula excerpt from VandyCRUISE
 * LesionToadSegmentation. It captures deterministic initial buffers and
 * statistics without pulling MIPAV UI, atlas, or evolution dependencies.
 */
public final class LesionSegmentationStateReferenceHarness {
    private static final int NX=3, NY=2, NZ=2, NC=2, K=3;
    private static final byte[] LABEL={1,2,3};
    private static final String[] TYPE={"obj","lesion","optimized_blackHole"};
    private static float[][][][] images() {
        float[][][][] a=new float[NC][NX][NY][NZ];
        for(int c=0;c<NC;c++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)
            a[c][x][y][z]=(float)(0.125+c*1.75+x*0.3+y*0.2+z*0.07);
        return a;
    }
    private static byte[][][] template() {
        byte[][][] a=new byte[NX][NY][NZ];
        int n=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)a[x][y][z]=(byte)((n++%3)==2?1:2);
        return a;
    }
    private static String floats(float[] a) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.isFinite(a[i])?Float.toString(a[i]):"null");}
        return s.append(']').toString();
    }
    private static String ints(int[] a) {
        StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();
    }
    private static String booleans(boolean[] a) {
        StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();
    }
    private static int signedBytePlusOne(int value) { return (byte)(((byte)value)+1); }
    public static void main(String[] args) {
        float[][][][] img=images();byte[][][] labels=template();
        byte[][][] seg=new byte[NX][NY][NZ],skeleton=new byte[NX][NY][NZ];
        short[][][] relation=new short[NX][NY][NZ];float[][][][] mem=new float[NX][NY][NZ][K];
        boolean[][][][] available=new boolean[NX][NY][NZ][K];float[] volume=new float[K];
        float[][] centroid=new float[NC][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            seg[x][y][z]=-1;skeleton[x][y][z]=-1;relation[x][y][z]=-1;
            for(int k=0;k<K;k++) if(labels[x][y][z]==LABEL[k]) {
                volume[k]++;seg[x][y][z]=(byte)k;skeleton[x][y][z]=(byte)k;relation[x][y][z]=(short)k;
                for(int c=0;c<NC;c++)centroid[c][k]+=img[c][x][y][z];
            }
        }
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++)centroid[c][k]/=volume[k];
        int n=NX*NY*NZ;
        float[] allImages=new float[NC*n];for(int c=0;c<NC;c++){int i=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)allImages[c*n+i++]=img[c][x][y][z];}
        int[] flatLabels=new int[n],flatSeg=new int[n],flatRelation=new int[n];int pos=0;
        float[] flatMem=new float[n*K];boolean[] flatAvailable=new boolean[n*K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            flatLabels[pos]=labels[x][y][z];flatSeg[pos]=seg[x][y][z];flatRelation[pos]=relation[x][y][z];pos++;
        }
        float[] flatCentroid=new float[NC*K];for(int c=0;c<NC;c++)for(int k=0;k<K;k++)flatCentroid[c*K+k]=centroid[c][k];
        float[] classSimilarity=new float[K*K],lesionSimilarity=new float[K],blackHoleSimilarity=new float[K];
        boolean hasLesion=false;int lesionClass=-1;
        for(int k=0;k<K;k++){
            String type=TYPE[k];
            if(type.equals("lesion")||type.equals("blackhole")||type.equals("optimized_lesion")||type.equals("optimized_blackHole"))hasLesion=true;
            if(type.equals("lesion")||type.equals("blackhole")||type.equals("optimized_lesion"))lesionClass=k;
            lesionSimilarity[k]=blackHoleSimilarity[k]=1.0f;
            for(int m=0;m<K;m++)classSimilarity[k*K+m]=k==m?0.0f:1.0f;
        }
        if(hasLesion&&lesionClass>=0)lesionSimilarity[lesionClass]=blackHoleSimilarity[lesionClass]=0.0f;
        // A second source-formula case confirms the defaults when the atlas has no lesion types.
        String[] noLesionTypes={"obj","mask","out"};
        boolean noHasLesion=false;int noLesionClass=-1;
        for(int k=0;k<K;k++)if(noLesionTypes[k].equals("lesion")||noLesionTypes[k].equals("blackhole")||noLesionTypes[k].equals("optimized_lesion")||noLesionTypes[k].equals("optimized_blackHole"))noHasLesion=true;
        for(int k=0;k<K;k++)if(noLesionTypes[k].equals("lesion")||noLesionTypes[k].equals("blackhole")||noLesionTypes[k].equals("optimized_lesion"))noLesionClass=k;
        float[] noLesionSimilarity=new float[K],noBlackHoleSimilarity=new float[K];
        java.util.Arrays.fill(noLesionSimilarity,1.0f);java.util.Arrays.fill(noBlackHoleSimilarity,1.0f);
        if(noHasLesion&&noLesionClass>=0)noLesionSimilarity[noLesionClass]=noBlackHoleSimilarity[noLesionClass]=0.0f;
        String json="{\"dimensions\":["+NX+","+NY+","+NZ+"],\"channels\":"+NC+",\"classes\":"+K
            +",\"labels\":"+ints(flatLabels)+",\"images\":"+floats(allImages)
            +",\"segmentation\":"+ints(flatSeg)+",\"skeleton\":"+ints(flatSeg)+",\"relation_map\":"+ints(flatRelation)
            +",\"memberships\":"+floats(flatMem)+",\"available\":"+booleans(flatAvailable)
            +",\"class_volume\":"+floats(volume)+",\"image_centroid\":"+floats(flatCentroid)
            +",\"absent_class_centroid_is_nan\":"+(!Float.isFinite(centroid[0][2]))
            +",\"class_similarity\":"+floats(classSimilarity)+",\"lesion_similarity\":"+floats(lesionSimilarity)
            +",\"black_hole_similarity\":"+floats(blackHoleSimilarity)
            +",\"no_lesion_similarity\":"+floats(noLesionSimilarity)+",\"no_lesion_black_hole_similarity\":"+floats(noBlackHoleSimilarity)
            +",\"ventricle_sentinel\":"+signedBytePlusOne(9)+",\"gm_sentinel\":"+signedBytePlusOne(5)
            +",\"brainstem_sentinel\":"+signedBytePlusOne(5)+",\"lesion_class\":"+lesionClass+",\"no_lesion_class\":"+noLesionClass+"}";
        System.out.println(json);
    }
}
