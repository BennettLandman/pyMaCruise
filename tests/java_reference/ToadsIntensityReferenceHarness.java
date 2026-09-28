/*
 * Test-only MaCRUISE TOADS intensity-normalization reference.
 *
 * The robustMinimum/robustMaximum routines below are copied from
 * VandyCRUISE/src/edu/jhmi/rad/medic/libraries/ImageFunctions.java
 * (author: Pierre-Louis Bazin, version July 2004). MultipleToadPreprocess
 * calls these with ratio 0.01 and scales 2, then divides each image value by
 * (Ihigh-Ilow) without subtracting Ilow. This harness has no MIPAV/JIST
 * dependency and emits reference values for a deterministic masked phantom.
 */
import java.util.Locale;

public final class ToadsIntensityReferenceHarness {
    private static final int NX = 13;
    private static final int NY = 11;
    private static final int NZ = 9;

    private static float[] phantom() {
        int size = NX * NY * NZ;
        float[] image = new float[size];
        for (int i = 0; i < size; i++) {
            int x = i % NX;
            int y = (i / NX) % NY;
            int z = i / (NX * NY);
            boolean inside = x >= 1 && x < NX - 1 && y >= 1 && y < NY - 1
                    && z >= 1 && z < NZ - 1;
            float value = 0.5f + ((i * 37) % 501) / 10.0f;
            if (i == 351 || i == 527 || i == 703) value = 300.0f + i;
            image[i] = inside ? value : 0.0f;
        }
        return image;
    }

    private static float robustMinimum(float[] image, float ratio, int scales) {
        float min = image[0], max = image[0];
        for (float value : image) { if (value > max) max = value; if (value < min) min = value; }
        ratio = ratio * image.length;
        float result = min;
        for (int t = 0; t < scales; t++) {
            float[] bins = histogram(image, min, max);
            float count = 0.0f;
            int n = 0;
            while (count < ratio && n < 10) count += bins[n++];
            result = min + (n - 0.5f) / 10.0f * (max - min);
            float oldMin = min;
            min = oldMin + (n - 1.0f) / 10.0f * (max - oldMin);
            max = oldMin + n / 10.0f * (max - oldMin);
            ratio = ratio - (count - bins[n - 1]);
        }
        return result;
    }

    private static float robustMaximum(float[] image, float ratio, int scales) {
        float min = image[0], max = image[0];
        for (float value : image) { if (value > max) max = value; if (value < min) min = value; }
        ratio = ratio * image.length;
        float result = max;
        for (int t = 0; t < scales; t++) {
            float[] bins = histogram(image, min, max);
            float count = 0.0f;
            int n = 10;
            while (count < ratio && n > 0) count += bins[--n];
            result = min + (n + 0.5f) / 10.0f * (max - min);
            float oldMin = min;
            min = oldMin + n / 10.0f * (max - oldMin);
            max = oldMin + (n + 1.0f) / 10.0f * (max - oldMin);
            ratio = ratio - (count - bins[n]);
        }
        return result;
    }

    private static float[] histogram(float[] image, float min, float max) {
        float[] bins = new float[10];
        for (float value : image) {
            if (value >= min && value <= min + 1.0f / 10.0f * (max - min)) bins[0]++;
            for (int n = 1; n < 10; n++) {
                if (value > min + n / 10.0f * (max - min)
                        && value <= min + (n + 1.0f) / 10.0f * (max - min)) bins[n]++;
            }
        }
        return bins;
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float[] image = phantom();
        float low = robustMinimum(image, 0.01f, 2);
        float high = robustMaximum(image, 0.01f, 2);
        float span = high - low;
        System.out.print("{\"shape\":[13,11,9],\"low\":" + low + ",\"high\":" + high + ",\"normalized\":[");
        for (int i = 0; i < image.length; i++) {
            if (i > 0) System.out.print(",");
            System.out.print(image[i] / span);
        }
        System.out.println("]}");
    }
}
