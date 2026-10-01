import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Small, independently compiled excerpt of MedicAlgorithmLesionToads atlas
 * selection. The switch follows the supplied VandyCRUISE source at lines
 * 405-422. Descriptor interpretation follows LesionToadSegmentation lines
 * 221-250: an atlas class with topology type "lesion" enables hasLesion.
 * This pins selection/configuration only; it does not execute LesionTOADS.
 */
public final class MaCRUISEplusLesionAtlasSelectionReferenceHarness {
    private static final class AtlasClass {
        final int slot;
        final String name;
        final int label;
        final String type;

        AtlasClass(int slot, String name, int label, String type) {
            this.slot = slot;
            this.name = name;
            this.label = label;
            this.type = type;
        }
    }

    private static Path selectAtlas(
            int atlasSelect,
            boolean mprage,
            boolean flair,
            Path withLesion,
            Path noLesionT1Flair,
            Path noLesionT1Only) {
        switch (atlasSelect) {
            case 0:
                return withLesion;
            case 1:
                if (mprage && flair) {
                    return noLesionT1Flair;
                } else {
                    return noLesionT1Only;
                }
            default:
                return withLesion;
        }
    }

    private static List<AtlasClass> parseClasses(Path descriptor) throws Exception {
        List<AtlasClass> classes = new ArrayList<>();
        boolean inClasses = false;
        for (String raw : Files.readAllLines(descriptor)) {
            String line = raw.trim();
            if (line.startsWith("Structures:")) {
                inClasses = true;
                continue;
            }
            if (line.startsWith("Topology Atlas:")) {
                break;
            }
            if (!inClasses || line.isEmpty()) {
                continue;
            }
            String[] fields = line.split("\\s+");
            if (fields.length >= 3 && fields[1].matches("-?\\d+")) {
                classes.add(new AtlasClass(
                        classes.size(), fields[0], Integer.parseInt(fields[1]), fields[2]));
            }
        }
        return classes;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 6) {
            throw new IllegalArgumentException(
                    "usage: <atlas-index> <mprage-present> <flair-present> "
                            + "<with-lesion> <no-lesion-t1-flair> <no-lesion-t1-only>");
        }
        Path selected = selectAtlas(
                Integer.parseInt(args[0]),
                Boolean.parseBoolean(args[1]),
                Boolean.parseBoolean(args[2]),
                Path.of(args[3]), Path.of(args[4]), Path.of(args[5]));
        String classEvidence = "not-inspected";
        if (Files.isRegularFile(selected)) {
            List<AtlasClass> classes = parseClasses(selected);
            boolean hasLesion = false;
            for (AtlasClass item : classes) {
                if (item.type.equals("lesion") || item.type.equals("blackhole")
                        || item.type.equals("optimized_lesion")
                        || item.type.equals("optimized_blackHole")) {
                    hasLesion = true;
                }
                if (item.name.equals("OuterVentricular-WM")) {
                    classEvidence = item.slot + "|" + item.name + "|" + item.label
                            + "|" + item.type + "|hasLesion=" + hasLesion;
                }
            }
        }
        System.out.println(selected.toAbsolutePath().normalize() + "\t" + classEvidence);
    }
}
