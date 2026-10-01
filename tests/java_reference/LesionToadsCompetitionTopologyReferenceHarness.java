import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Random;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;

/**
 * Selected LesionToadSegmentation.isHomeomorphicSegmentation excerpt at the
 * VandyCRUISE wrapper's explicit 18/6 connectivity.
 */
public final class LesionToadsCompetitionTopologyReferenceHarness {
    private static final CriticalPointLUT LUT =
            new CriticalPointLUT("/edu/jhmi/rad/medic/structures/critical186LUT.raw.gz", 200);

    private static boolean isHomeomorphic(short[][][] segmentation, int x, int y,
            int z, short candidate) {
        if (segmentation[x][y][z] == candidate) return true;
        boolean[][][] obj = new boolean[3][3][3];

        // Directly preserve the selected LesionToadSegmentation method body.
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
            obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == candidate;
        obj[1][1][1] = true;
        if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;

        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
            obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == segmentation[x][y][z];
        obj[1][1][1] = false;
        if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;

        int count = 0;
        short[] labels = new short[26];
        for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
            if ((i*i+j*j+l*l>0) && segmentation[x+i][y+j][z+l] != candidate
                    && segmentation[x+i][y+j][z+l] != segmentation[x][y][z]) {
                boolean found = false;
                for (int n=0;n<count;n++) if (segmentation[x+i][y+j][z+l] == labels[n]) {
                    found = true; break;
                }
                if (!found) labels[count++] = segmentation[x+i][y+j][z+l];
            }
        }

        for (int n=0;n<count;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
                obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == segmentation[x][y][z]
                        || segmentation[x+i][y+j][z+l] == labels[n];
            obj[1][1][1] = false;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
                obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == candidate
                        || segmentation[x+i][y+j][z+l] == labels[n];
            obj[1][1][1] = true;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) for (int m=n+1;m<count;m++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                short label = segmentation[x+i][y+j][z+l];
                obj[1+i][1+j][1+l] = label == segmentation[x][y][z]
                        || label == labels[n] || label == labels[m];
            }
            obj[1][1][1] = false;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) for (int m=n+1;m<count;m++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                short label = segmentation[x+i][y+j][z+l];
                obj[1+i][1+j][1+l] = label == candidate
                        || label == labels[n] || label == labels[m];
            }
            obj[1][1][1] = true;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }
        return true;
    }

    private static String encode(short[][][] neighborhood) {
        StringBuilder out = new StringBuilder(27);
        for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
            out.append((char)('0' + neighborhood[x][y][z]));
        return out.toString();
    }

    public static void main(String[] args) {
        PrintStream stdout = System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));
        boolean loaded = LUT.loadCompressedPattern();
        System.setOut(stdout);
        if (!loaded) throw new IllegalStateException("could not load 18/6 LUT");
        Random random = new Random(18618);
        System.out.print("{\"source_method\":\"LesionToadSegmentation.isHomeomorphicSegmentation\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"source_method_lines\":\"1611-1725\",\"connectivity\":\"18/6\",\"lut_source_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"lut_asset_sha256\":\"cbf1e339d78a0a2c673c939ee2c8a30ee3a35bfd9603d07ed2d48864bb0e88be\",\"cases\":[");
        for (int caseIndex=0;caseIndex<140;caseIndex++) {
            short[][][] neighborhood = new short[3][3][3];
            for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
                neighborhood[x][y][z] = (short)random.nextInt(5);
            short candidate = (short)random.nextInt(5);
            if (caseIndex == 0) candidate = neighborhood[1][1][1];
            short[][][] volume = new short[5][5][5];
            for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
                volume[x+1][y+1][z+1] = neighborhood[x][y][z];
            boolean allowed = isHomeomorphic(volume,2,2,2,candidate);
            if (caseIndex>0) System.out.print(",");
            System.out.printf("{\"neighborhood\":\"%s\",\"candidate\":%d,\"allowed\":%s}",
                    encode(neighborhood), candidate, allowed ? "true" : "false");
        }
        System.out.println("]}");
    }
}
