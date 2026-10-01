/**
 * Wrapper-order reference for field attachment before initial TOADS
 * computeExactMemberships. Numerical membership bits are pinned separately by
 * ToadsDefaultAtlasMembershipReferenceHarness.
 */
public final class ToadsInitialFieldedMembershipReferenceHarness {
    private static void addInhomogeneityCorrection(float[][][][] field,
            float[][][] field_, int c, boolean[] useField) {
        field[c] = field_;
        useField[0] = true;
    }

    public static void main(String[] args) {
        float[][][][] fields = new float[2][][][];
        boolean[] useField = {false};
        float[][][] t1 = {{{1.0f}}};
        float[][][] flair = {{{1.0f}}};
        addInhomogeneityCorrection(fields, t1, 0, useField);
        addInhomogeneityCorrection(fields, flair, 1, useField);
        String[] order = {
            "addInhomogeneityCorrection[channel=0]",
            "addInhomogeneityCorrection[channel=1]",
            "computeRelations",
            "computeBoundaryDistances",
            "computeClosestLabels",
            "computeVariances",
            "computeExactMemberships"
        };
        boolean references = fields[0] == t1 && fields[1] == flair;
        System.out.print("{\"use_field\":" + useField[0]
                + ",\"channel_references_preserved\":" + references
                + ",\"stage_order\":[");
        for (int i = 0; i < order.length; i++) {
            if (i > 0) System.out.print(',');
            System.out.print("\"" + order[i] + "\"");
        }
        System.out.println("]}");
    }
}
