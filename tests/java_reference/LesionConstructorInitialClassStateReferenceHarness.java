/**
 * Pinned initialization expressions from
 * VandyCRUISE MultipleToadSegmentation constructor, initial values loop.
 * Source: toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java
 * The harness uses a deterministic transformed-template byte volume.
 */
public final class LesionConstructorInitialClassStateReferenceHarness {
    private static String ints(byte[][][] volume) {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int x=0;x<volume.length;x++) for (int y=0;y<volume[x].length;y++) for (int z=0;z<volume[x][y].length;z++) {
            if (!first) out.append(','); first=false; out.append(volume[x][y][z]);
        }
        return out.append(']').toString();
    }
    private static String floats(float[][][] volume) {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int x=0;x<volume.length;x++) for (int y=0;y<volume[x].length;y++) for (int z=0;z<volume[x][y].length;z++) {
            if (!first) out.append(','); first=false;
            out.append(Integer.toUnsignedLong(Float.floatToIntBits(volume[x][y][z])));
        }
        return out.append(']').toString();
    }
    private static String bools(boolean[][][][] volume) {
        StringBuilder out = new StringBuilder("[");
        boolean first = true;
        for (int x=0;x<volume.length;x++) for (int y=0;y<volume[x].length;y++) for (int z=0;z<volume[x][y].length;z++)
            for (int k=0;k<volume[x][y][z].length;k++) {
                if (!first) out.append(','); first=false; out.append(volume[x][y][z][k]);
            }
        return out.append(']').toString();
    }
    public static void main(String[] args) {
        final int nx=4, ny=4, nz=4, classes=3;
        final byte EMPTY=-1;
        final byte[] templateLabel={1,5,24};
        final byte[][][] classification=new byte[nx][ny][nz];
        final byte[][][] segmentation=new byte[nx][ny][nz];
        final byte[][][] skeleton=new byte[nx][ny][nz];
        final byte[][][] relationMap=new byte[nx][ny][nz];
        final boolean[][][][] available=new boolean[nx][ny][nz][classes];
        final float[][][] ordering=new float[nx][ny][nz];
        final float[][][] lesionMembership=new float[nx][ny][nz];
        final float[][][] outlierMembership=new float[nx][ny][nz];
        final float[] volume=new float[classes];
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++)
            classification[x][y][z]=templateLabel[(x+2*y+z)%classes];
        for (int k=0;k<classes;k++) volume[k]=0;
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) {
            ordering[x][y][z]=0.0f;
            segmentation[x][y][z]=EMPTY;
            skeleton[x][y][z]=EMPTY;
            relationMap[x][y][z]=EMPTY;
            for(byte k=0;k<classes;k++) {
                if(classification[x][y][z]==templateLabel[k]) {
                    volume[k]++;
                    segmentation[x][y][z]=k;
                    skeleton[x][y][z]=k;
                    relationMap[x][y][z]=k;
                }
                available[x][y][z][k]=false;
            }
            lesionMembership[x][y][z]=0.25f;
            outlierMembership[x][y][z]=0.5f;
        }
        StringBuilder volumes=new StringBuilder("[");
        for(int k=0;k<classes;k++){if(k>0)volumes.append(',');volumes.append(Integer.toUnsignedLong(Float.floatToIntBits(volume[k])));}
        volumes.append(']');
        System.out.println("{\"dimensions\":[4,4,4],\"classes\":3,\"segmentation\":"+ints(segmentation)
            +",\"skeleton\":"+ints(skeleton)+",\"relation_map\":"+ints(relationMap)
            +",\"available\":"+bools(available)+",\"ordering_bits\":"+floats(ordering)
            +",\"lesion_membership_bits\":"+floats(lesionMembership)
            +",\"outlier_membership_bits\":"+floats(outlierMembership)+",\"class_volume_bits\":"+volumes+"}");
    }
}
