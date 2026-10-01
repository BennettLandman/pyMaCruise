/*
 * Test-only LesionTOADS preprocessing reference.
 *
 * robustMinimum/robustMaximum are copied from
 * toads-cruise/src/edu/jhmi/rad/medic/libraries/ImageFunctions.java.
 * The selected wrapper pins robustMaxMin=false. The clamp and ordinary
 * min/max follow the LesionToadPreprocess constructor; normalization follows
 * normalizeImages; robust bounds and low+prior*(high-low) centroid expressions
 * follow initialLesionCentroid and initialBlackHoleCentroid, on normalized
 * whole-grid inputs (the outer wrapper crop is not part of this harness).
 * This standalone harness uses float arithmetic and does not require MIPAV/JIST.
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
        float[] imageMin = new float[CHANNELS];
        float[] imageMax = new float[CHANNELS];
        float[] lesion = new float[CHANNELS];
        float[] blackHole = new float[CHANNELS];
        float[] croppedLow = new float[CHANNELS];
        float[] croppedHigh = new float[CHANNELS];
        float[] croppedLesion = new float[CHANNELS];
        float[] croppedBlackHole = new float[CHANNELS];
        for (int c = 0; c < CHANNELS; c++) {
            imageMin[c] = images[c][0];
            imageMax[c] = images[c][0];
            for (int i = 0; i < images[c].length; i++) {
                if (images[c][i] < 0.0f) {
                    images[c][i] = 0.0f;
                    imageMin[c] = 0.0f;
                } else {
                    if (images[c][i] > imageMax[c]) imageMax[c] = images[c][i];
                    if (images[c][i] < imageMin[c]) imageMin[c] = images[c][i];
                }
            }
            float span = imageMax[c] - imageMin[c];
            for (int i = 0; i < images[c].length; i++) {
                images[c][i] = (images[c][i] - imageMin[c]) / span;
            }
            low[c] = robustMinimum(images[c], 0.01f, 2);
            high[c] = robustMaximum(images[c], 0.01f, 2);
            float centroidSpan = high[c] - low[c];
            lesion[c] = low[c] + LESION_PRIORS[c] * centroidSpan;
            blackHole[c] = low[c] + BLACK_HOLE_PRIORS[c] * centroidSpan;

            // Model the selected wrapper's already-cropped, normalized image
            // at the centroid-method boundary. Cropping itself is not under
            // test here; use a fixed asymmetric interior region so robust
            // centroid bounds cannot accidentally reuse full-image bounds.
            float[] croppedNormalized = new float[(NX - 2) * (NY - 4) * (NZ - 2)];
            int index = 0;
            for (int z = 1; z < NZ - 1; z++) {
                for (int y = 2; y < NY - 2; y++) {
                    for (int x = 1; x < NX - 1; x++) {
                        croppedNormalized[index++] = images[c][x + NX * y + NX * NY * z];
                    }
                }
            }
            croppedLow[c] = robustMinimum(croppedNormalized, 0.01f, 2);
            croppedHigh[c] = robustMaximum(croppedNormalized, 0.01f, 2);
            float croppedSpan = croppedHigh[c] - croppedLow[c];
            croppedLesion[c] = croppedLow[c] + LESION_PRIORS[c] * croppedSpan;
            croppedBlackHole[c] = croppedLow[c] + BLACK_HOLE_PRIORS[c] * croppedSpan;
        }

        System.out.print("{\"source_sha256\":{\"LesionToadPreprocess.java\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",\"AlgorithmLesionToads.java\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"ImageFunctions.java\":\"be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892\"},\"call_order\":[\"LesionToadPreprocess.constructor\",\"AlgorithmLesionToads.normalizeImages\",\"findCroppingBoundaries\",\"cropImages\",\"initialCentroids\",\"initialLesionCentroid\",\"initialBlackHoleCentroid\"],\"reference_scope\":\"No-crop source excerpt: normalization then both auxiliary centroid methods over the full normalized grid\",\"robustMaxMin\":false,\"shape\":[" + NX + "," + NY + "," + NZ + "],\"channels\":[");
        for (int c = 0; c < CHANNELS; c++) {
            if (c > 0) System.out.print(',');
            System.out.print(array(images[c]));
        }
        System.out.print("],\"intensity_minima\":" + array(imageMin));
        System.out.print(",\"intensity_maxima\":" + array(imageMax));
        System.out.print(",\"centroid_robust_minima\":" + array(low));
        System.out.print(",\"centroid_robust_maxima\":" + array(high));
        System.out.print(",\"lesion_priors\":" + array(LESION_PRIORS));
        System.out.print(",\"black_hole_priors\":" + array(BLACK_HOLE_PRIORS));
        System.out.print(",\"lesion_centroid\":" + array(lesion));
        System.out.print(",\"black_hole_centroid\":" + array(blackHole));
        System.out.print(",\"cropped_centroid_robust_minima\":" + array(croppedLow));
        System.out.print(",\"cropped_centroid_robust_maxima\":" + array(croppedHigh));
        System.out.print(",\"cropped_lesion_centroid\":" + array(croppedLesion));
        System.out.println(",\"cropped_black_hole_centroid\":" + array(croppedBlackHole) + "}");
    }
}
