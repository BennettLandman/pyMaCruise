import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Random;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;

/** Calls the selected VandyCRUISE LUT key and lookup methods directly. */
public final class ToadsCriticalPointLutReferenceHarness {
    private static final int RANDOM_CASES = 64;

    private static String encode(boolean[][][] pattern) {
        StringBuilder bits = new StringBuilder(27);
        for (int x = 0; x < 3; x++) {
            for (int y = 0; y < 3; y++) {
                for (int z = 0; z < 3; z++) {
                    bits.append(pattern[x][y][z] ? '1' : '0');
                }
            }
        }
        return bits.toString();
    }

    private static void emit(CriticalPointLUT lut, boolean[][][] pattern,
            boolean first) {
        if (!first) System.out.print(",");
        int key = lut.keyFromPattern(pattern, 1, 1, 1);
        System.out.printf("{\"neighborhood\":\"%s\",\"key\":%d,\"regular\":%s}",
                encode(pattern), key, lut.get(key) ? "true" : "false");
    }

    public static void main(String[] args) {
        CriticalPointLUT lut = new CriticalPointLUT("/critical186LUT.raw.gz", 200);
        PrintStream output = System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));
        boolean loaded = lut.loadCompressedPattern();
        System.setOut(output);
        if (!loaded) {
            throw new IllegalStateException("Could not load bundled reference LUT");
        }
        System.out.print("{\"source\":\"VandyCRUISE CriticalPointLUT\",\"cases\":[");

        boolean[][][] pattern = new boolean[3][3][3];
        emit(lut, pattern, true);

        pattern[0][0][0] = true;
        emit(lut, pattern, false);

        pattern = new boolean[3][3][3];
        pattern[0][1][1] = true;
        pattern[2][1][1] = true;
        emit(lut, pattern, false);

        pattern = new boolean[3][3][3];
        for (int x = 0; x < 3; x++) {
            for (int y = 0; y < 3; y++) {
                for (int z = 0; z < 3; z++) pattern[x][y][z] = true;
            }
        }
        emit(lut, pattern, false);

        Random random = new Random(186);
        for (int n = 0; n < RANDOM_CASES; n++) {
            pattern = new boolean[3][3][3];
            for (int x = 0; x < 3; x++) {
                for (int y = 0; y < 3; y++) {
                    for (int z = 0; z < 3; z++) pattern[x][y][z] = random.nextBoolean();
                }
            }
            emit(lut, pattern, false);
        }
        System.out.println("]}");
    }
}
