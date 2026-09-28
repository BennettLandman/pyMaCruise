/*
 * Test-only reference for the selected TOADS descriptor intensity-prior
 * lookup. The alias mapping follows
 * toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java
 * method modalityId; the profile loop follows its Intensity Atlas parser.
 * The original atlas/parser source is by Pierre-Louis Bazin. This harness
 * uses only the descriptor text and Java float parsing; no MIPAV/JIST runtime.
 */
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class ToadsIntensityPriorReferenceHarness {
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
        if (args.length < 2) throw new IllegalArgumentException("descriptor modality...");
        List<String> lines = Files.readAllLines(Path.of(args[0]));
        int classes = -1;
        int intensityStart = -1;
        for (int i = 0; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.startsWith("Structures:")) classes = Integer.parseInt(line.split("\\s+")[1]);
            if (line.startsWith("Intensity Atlas:")) intensityStart = i + 1;
        }
        if (classes < 1 || intensityStart < 0) throw new IllegalArgumentException("missing sections");

        Map<String, float[]> profiles = new LinkedHashMap<>();
        for (int i = intensityStart; i < lines.size(); i++) {
            String[] fields = lines.get(i).trim().split("\\s+");
            if (fields.length == 1 && fields[0].isEmpty()) continue;
            String id = modalityId(fields[0]);
            if (id == null) break;
            if (fields.length != classes + 1) throw new IllegalArgumentException("class count mismatch");
            float[] values = new float[classes];
            for (int k = 0; k < classes; k++) values[k] = Float.parseFloat(fields[k + 1]);
            profiles.put(id, values);
        }

        List<String> outputModalities = new ArrayList<>();
        List<float[]> outputPriors = new ArrayList<>();
        for (int i = 1; i < args.length; i++) {
            String id = modalityId(args[i]);
            if (id == null || !profiles.containsKey(id)) throw new IllegalArgumentException("missing requested modality");
            outputModalities.add(args[i]);
            outputPriors.add(profiles.get(id));
        }
        StringBuilder output = new StringBuilder("{\"class_count\":" + classes + ",\"modalities\":[");
        for (int i = 0; i < outputModalities.size(); i++) {
            if (i > 0) output.append(',');
            output.append(escape(outputModalities.get(i)));
        }
        output.append("],\"priors\":[");
        for (int i = 0; i < outputPriors.size(); i++) {
            if (i > 0) output.append(',');
            output.append('[');
            for (int k = 0; k < outputPriors.get(i).length; k++) {
                if (k > 0) output.append(',');
                output.append(Float.toString(outputPriors.get(i)[k]));
            }
            output.append(']');
        }
        System.out.println(output.append("]}").toString());
    }
}
