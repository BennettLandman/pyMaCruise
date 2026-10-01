/**
 * Narrow selected-source reference for MultipleToadSegmentation's field setter.
 * It intentionally tests assignment/state only, not field fitting or consumers.
 */
public final class ToadsInhomogeneityFieldAttachmentReferenceHarness {
    private static void addInhomogeneityCorrection(float[][][][] field,
            float[][][] field_, int c, boolean[] useField) {
        field[c] = field_;
        useField[0] = true;
    }

    public static void main(String[] args) {
        float[][][][] fields = new float[2][][][];
        boolean[] useField = {false};
        float[][][] first = {{{1.0f, 1.25f}}, {{0.75f, 1.5f}}};
        float[][][] second = {{{2.0f, 2.25f}}, {{1.75f, 2.5f}}};
        boolean initialFalse = !useField[0] && fields[0] == null && fields[1] == null;
        addInhomogeneityCorrection(fields, first, 1, useField);
        boolean assignedByIdentity = fields[1] == first;
        boolean otherChannelUntouched = fields[0] == null;
        boolean enabledAfterAssignment = useField[0];
        first[0][0][0] = 3.0f;
        boolean mutationVisible = fields[1][0][0][0] == 3.0f;
        addInhomogeneityCorrection(fields, second, 1, useField);
        boolean replacementByIdentity = fields[1] == second;
        addInhomogeneityCorrection(fields, null, 0, useField);
        boolean nullAcceptedAndEnabled = fields[0] == null && useField[0];
        System.out.println("{\"initial_false\":" + initialFalse
                + ",\"assigned_by_identity\":" + assignedByIdentity
                + ",\"other_channel_untouched\":" + otherChannelUntouched
                + ",\"enabled_after_assignment\":" + enabledAfterAssignment
                + ",\"mutation_visible\":" + mutationVisible
                + ",\"replacement_by_identity\":" + replacementByIdentity
                + ",\"null_accepted_and_enabled\":" + nullAcceptedAndEnabled + "}");
    }
}
