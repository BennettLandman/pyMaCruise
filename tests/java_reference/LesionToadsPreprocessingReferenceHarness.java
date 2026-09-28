/*
 * Test-only LesionTOADS preprocessing reference.
 *
 * robustMinimum/robustMaximum are copied from
 * toads-cruise/src/edu/jhmi/rad/medic/libraries/ImageFunctions.java.
 * The clamp, normalization, and lesion/black-hole centroid expressions follow
 * LesionToadPreprocess.java (constructor, normalizeImages,
 * initialLesionCentroid, and initialBlackHoleCentroid). This standalone
 * harness uses float arithmetic and does not require MIPAV or JIST.
 */
import java.util.Locale;

public final class LesionToadsPreprocessingReferenceHarness {
    private static final int NX = 9;
    private static final int NY = 8;
    private static final int NZ = 7;
    private static final int CHANNELS = 2;
    private static final float[] LESION_PRIORS = {0.4f, 1.2f};
    private static final float[] BLACK_HOLE_PRIORS = {0.1f, 0.6f};

    private static float[][] phantom() {
        int size = NX * NY * NZ;
        float[][] images = new float[CHANNELS][size];
        for (int i = 0; i < size; i++) {
            images[0][i] = 4.25f + ((i * 19 + 3) % 137) / 16.0f;
            images[1][i] = 11.5f + ((i * 29 + 5) % 211) / 12.0f;
        }
        for (int i = 0; i < size; i += 61) {
            images[0][i] = -5.0f - i;
            images[1][i] = -9.0f - i;
        }
        images[0][37] = 600.0f;
        images[0][414] = 1200.0f;
        images[1][73] = 900.0f;
        images[1][422] = 1600.0f;
        return images;
    }

    private static float robustMinimum(float[] image, float ratio, int scales) {
        float min = image[0], max = image[0];
        for (float value : image) {
            if (value > max) max = value;
            if (value < min) min = value;
        }
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
        for (float value : image) {
            if (value > max) max = value;
            if (value < min) min = value;
        }
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

    private static String array(float[] values) {
        StringBuilder text = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) text.append(',');
            text.append(Float.toString(values[i]));
        }
        return text.append(']').toString();
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        float[][] images = phantom();
        float[] low = new float[CHANNELS];
        float[] high = new float[CHANNELS];
        float[] lesion = new float[CHANNELS];
        float[] blackHole = new float[CHANNELS];
        for (int c = 0; c < CHANNELS; c++) {
            for (int i = 0; i < images[c].length; i++) {
                if (images[c][i] < 0.0f) images[c][i] = 0.0f;
            }
            low[c] = robustMinimum(images[c], 0.01f, 2);
            high[c] = robustMaximum(images[c], 0.01f, 2);
            float span = high[c] - low[c];
            float lowOverSpan = low[c] / span;
            lesion[c] = lowOverSpan + LESION_PRIORS[c];
            blackHole[c] = lowOverSpan + BLACK_HOLE_PRIORS[c];
            for (int i = 0; i < images[c].length; i++) images[c][i] = images[c][i] / span;
        }

        System.out.print("{\"shape\":[" + NX + "," + NY + "," + NZ + "],\"channels\":[");
        for (int c = 0; c < CHANNELS; c++) {
            if (c > 0) System.out.print(',');
            System.out.print(array(images[c]));
        }
        System.out.print("],\"robust_minima\":" + array(low));
        System.out.print(",\"robust_maxima\":" + array(high));
        System.out.print(",\"lesion_priors\":" + array(LESION_PRIORS));
        System.out.print(",\"black_hole_priors\":" + array(BLACK_HOLE_PRIORS));
        System.out.print(",\"lesion_centroid\":" + array(lesion));
        System.out.println(",\"black_hole_centroid\":" + array(blackHole) + "}");
    }
}
