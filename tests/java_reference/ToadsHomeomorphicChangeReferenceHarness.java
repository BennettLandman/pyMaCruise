import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Random;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;

/**
 * Executes the selected isHomeomorphicSegmentation method body on synthetic
 * local neighborhoods, with checkComposed=false as configured by the source.
 */
public final class ToadsHomeomorphicChangeReferenceHarness {
    private static final CriticalPointLUT LUT = new CriticalPointLUT("/critical186LUT.raw.gz", 200);
    private static final int RANDOM_CASES = 128;

    private static boolean isHomeomorphic(short[][][] segmentation, short x, short y,
            short z, short k) {
        if (segmentation[x][y][z] == k) return true;
        boolean[][][] obj = new boolean[3][3][3];
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
            obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == k;
        }
        obj[1][1][1] = true;
        if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;

        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
            obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == segmentation[x][y][z];
        }
        obj[1][1][1] = false;
        if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;

        int nConfiguration = 0;
        short[] labels = new short[26];
        for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
            if ((i*i+j*j+l*l>0)
                    && segmentation[x+i][y+j][z+l] != k
                    && segmentation[x+i][y+j][z+l] != segmentation[x][y][z]) {
                boolean found = false;
                for (int n=0;n<nConfiguration;n++) {
                    if (segmentation[x+i][y+j][z+l] == labels[n]) { found = true; break; }
                }
                if (!found) labels[nConfiguration++] = segmentation[x+i][y+j][z+l];
            }
        }

        // Pair relations: old object, then new object.
        for (int n=0;n<nConfiguration;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == segmentation[x][y][z]
                        || segmentation[x+i][y+j][z+l] == labels[n];
            }
            obj[1][1][1] = false;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<nConfiguration;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                obj[1+i][1+j][1+l] = segmentation[x+i][y+j][z+l] == k
                        || segmentation[x+i][y+j][z+l] == labels[n];
            }
            obj[1][1][1] = true;
            if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
        }

        // Triplet relations: old object, then new object.
        for (int n=0;n<nConfiguration;n++) {
            for (int m=n+1;m<nConfiguration;m++) {
                for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                    short label = segmentation[x+i][y+j][z+l];
                    obj[1+i][1+j][1+l] = label == segmentation[x][y][z]
                            || label == labels[n] || label == labels[m];
                }
                obj[1][1][1] = false;
                if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
            }
        }
        for (int n=0;n<nConfiguration;n++) {
            for (int m=n+1;m<nConfiguration;m++) {
                for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                    short label = segmentation[x+i][y+j][z+l];
                    obj[1+i][1+j][1+l] = label == k || label == labels[n] || label == labels[m];
                }
                obj[1][1][1] = true;
                if (!LUT.get(LUT.keyFromPattern(obj,1,1,1))) return false;
            }
        }
        return true;
    }

    private static String encode(short[][][] neighborhood) {
        StringBuilder result = new StringBuilder(27);
        for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
            result.append((char)('0' + neighborhood[x][y][z]));
        return result.toString();
    }

    private static void emit(short[][][] neighborhood, short candidate, boolean first) {
        short[][][] segmentation = new short[5][5][5];
        for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
            segmentation[x+1][y+1][z+1] = neighborhood[x][y][z];
        boolean allowed = isHomeomorphic(segmentation, (short)2, (short)2, (short)2, candidate);
        if (!first) System.out.print(",");
        System.out.printf("{\"neighborhood\":\"%s\",\"candidate\":%d,\"allowed\":%s}",
                encode(neighborhood), candidate, allowed ? "true" : "false");
    }

    public static void main(String[] args) {
        PrintStream output = System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));
        if (!LUT.loadCompressedPattern()) throw new IllegalStateException("Could not load 18/6 LUT");
        System.setOut(output);
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.isHomeomorphicSegmentation excerpt\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"source_lines\":\"1554-1673\",\"lut_source_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"check_composed\":false,\"cases\":[");

        short[][][] neighborhood = new short[3][3][3];
        emit(neighborhood, (short)1, true);
        neighborhood[1][1][1] = 1;
        emit(neighborhood, (short)1, false);

        Random random = new Random(1866);
        for (int n=0;n<RANDOM_CASES;n++) {
            neighborhood = new short[3][3][3];
            for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++)
                neighborhood[x][y][z] = (short)random.nextInt(4);
            emit(neighborhood, (short)random.nextInt(4), false);
        }
        System.out.println("]}");
    }
}
