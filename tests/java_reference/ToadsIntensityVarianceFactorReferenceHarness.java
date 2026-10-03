/*
 * Standalone reference harness for
 * MultipleToadDeformableAtlas.getIntensityVariancePriors(String[], int).
 * The selected source initializes every variance entry to 1.0f at lines
 * 939-942, parses only the "Intensity Variance Factor" section at lines
 * 1114-1124, and returns parsed rows for known modalities while leaving
 * unknown modalities at the initialized ones (lines 270-279). modalityId is
 * mirrored below. This isolates that source behavior without MIPAV/JIST.
 */
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class ToadsIntensityVarianceFactorReferenceHarness {
    private static String modalityId(String name) {
        switch (name) {
            case "T1": case "T1_SPGR": return "T1_SPGR";
            case "T2": return "T2";
            case "FLAIR": return "FLAIR";
            case "MPRAGE": case "T1_MPRAGE": return "T1_MPRAGE";
            case "MPRAGE_RAW": case "T1_RAW": return "T1_RAW";
            case "PD": return "PD";
            case "PDFSE": return "PDFSE";
            case "DIR": return "DIR";
            default: return null;
        }
    }

    private static String escape(String value) {
        return "\"" + value.replace("\\", "\\\\").replace("\"", "\\\"") + "\"";
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 1) throw new IllegalArgumentException("descriptor [modality...]");
        List<String> lines = Files.readAllLines(Path.of(args[0]));
        int classes = -1;
        int varianceStart = -1;
        for (int i = 0; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.startsWith("Structures:")) classes = Integer.parseInt(line.split("\\s+")[1]);
            if (line.startsWith("Intensity Variance Factor:")) varianceStart = i + 1;
        }
        if (classes < 1) throw new IllegalArgumentException("missing or invalid Structures section");

        Map<String, float[]> profiles = new LinkedHashMap<>();
        if (varianceStart >= 0) {
            for (int i = varianceStart; i < lines.size(); i++) {
                String line = lines.get(i).trim();
                if (line.isEmpty()) continue;
                String[] fields = line.split("\\s+");
                String id = modalityId(fields[0]);
                if (id == null) break;
                if (fields.length != classes + 1) throw new IllegalArgumentException("class count mismatch");
                float[] values = new float[classes];
                for (int k = 0; k < classes; k++) values[k] = Float.parseFloat(fields[k + 1]);
                profiles.put(id, values);
            }
        }

        List<String> outputModalities = new ArrayList<>();
        List<float[]> output = new ArrayList<>();
        for (int i = 1; i < args.length; i++) {
            outputModalities.add(args[i]);
            float[] row = new float[classes];
            java.util.Arrays.fill(row, 1.0f);
            String id = modalityId(args[i]);
            if (id != null && profiles.containsKey(id)) row = profiles.get(id);
            output.add(row);
        }
        StringBuilder json = new StringBuilder("{\"class_count\":" + classes + ",\"modalities\":[");
        for (int i = 0; i < outputModalities.size(); i++) {
            if (i > 0) json.append(',');
            json.append(escape(outputModalities.get(i)));
        }
        json.append("],\"factors\":[");
        for (int i = 0; i < output.size(); i++) {
            if (i > 0) json.append(',');
            json.append('[');
            for (int k = 0; k < classes; k++) {
                if (k > 0) json.append(',');
                json.append(Float.toString(output.get(i)[k]));
            }
            json.append(']');
        }
        System.out.println(json.append("]}").toString());
    }
}
