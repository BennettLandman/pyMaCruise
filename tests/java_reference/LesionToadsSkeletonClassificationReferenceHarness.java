/** Selected-method excerpt of VandyCRUISE LesionToadSegmentation.updateClassificationFromSkeleton. */
public final class LesionToadsSkeletonClassificationReferenceHarness {
    private static int updateClassificationFromSkeleton(byte[][][] classification,
            byte[] templateLabel, byte[][][] skeleton) {
        int nx = classification.length, ny = classification[0].length, nz = classification[0][0].length;
        int classDiff = 0;
        byte prev;
        for (int x = 0; x < nx; x++) for (int y = 0; y < ny; y++) for (int z = 0; z < nz; z++) {
            prev = classification[x][y][z];
            classification[x][y][z] = 0;
            if (skeleton[x][y][z] != (byte)-1)
                classification[x][y][z] = templateLabel[skeleton[x][y][z]];
            if (prev != classification[x][y][z]) classDiff++;
        }
        return classDiff;
    }

    public static void main(String[] args) {
        byte[] templateLabel = {(byte)-10, (byte)20, (byte)-128, (byte)99};
        byte[] initial = {(byte)5, (byte)-10, (byte)20, (byte)0, (byte)99,
            (byte)0, (byte)2, (byte)-10, (byte)-7, (byte)20};
        byte[] skeletonValues = {(byte)-1, (byte)0, (byte)1, (byte)2, (byte)3,
            (byte)-1, (byte)2, (byte)0, (byte)3, (byte)1};
        byte[][][] classification = new byte[10][1][1];
        byte[][][] skeleton = new byte[10][1][1];
        for (int x = 0; x < 10; x++) {
            classification[x][0][0] = initial[x];
            skeleton[x][0][0] = skeletonValues[x];
        }
        int changed = updateClassificationFromSkeleton(classification, templateLabel, skeleton);
        System.out.print("{\"source_method\":\"LesionToadSegmentation.updateClassificationFromSkeleton\",");
        System.out.print("\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",");
        System.out.print("\"template_labels\":[-10,20,-128,99],\"initial\":[");
        for (int i = 0; i < 10; i++) { if (i > 0) System.out.print(","); System.out.print(initial[i]); }
        System.out.print("],\"skeleton\":[");
        for (int i = 0; i < 10; i++) { if (i > 0) System.out.print(","); System.out.print(skeletonValues[i]); }
        System.out.print("],\"classification\":[");
        for (int x = 0; x < 10; x++) { if (x > 0) System.out.print(","); System.out.print(classification[x][0][0]); }
        System.out.println("],\"changed\":" + changed + "}");
    }
}
