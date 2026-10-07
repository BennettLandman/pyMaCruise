import java.io.File;

import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.tgdm.ProfileTGDM;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.vanderbilt.masi.plugins.CRUISE.tgdm.MaGenericTGDM;

/** Isolate the selected Java outer-surface call using its Java central field. */
public final class TGDMOuterStageReferenceHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 3) {
            throw new IllegalArgumentException("usage: TGDMOuterStageReferenceHarness INPUT_DIR CORE_REFERENCE_DIR OUTPUT_DIR");
        }
        File input = new File(args[0]);
        File core = new File(args[1]);
        File output = new File(args[2]);
        if (!output.isDirectory() && !output.mkdirs())
            throw new java.io.IOException("cannot create output directory: " + output);
        ImageDataFloat central = new ImageDataFloat(
                TGDMCoreReferenceHarness.readNifti(new File(core, "central.nii.gz")));
        ImageDataFloat gray = new ImageDataFloat(
                TGDMCoreReferenceHarness.readNifti(new File(input, "gray_matter.nii.gz")));
        ImageDataFloat skeleton = new ImageDataFloat(
                TGDMCoreReferenceHarness.readNifti(new File(input, "ace_skeleton.nii.gz")));
        CapturingTGDM tgdm = new CapturingTGDM();
        ImageDataFloat result = tgdm
                .solveOuterSurface(central, null, gray, null, skeleton,
                        1.0f, 0.0f, 1.0f, 0.4f, 1);
        TGDMCoreReferenceHarness.writeNifti(
                new File(output, "outer_from_java_central.nii.gz"), result.toArray3d());
        TGDMCoreReferenceHarness.writeNifti(
                new File(output, "outer_before_reinit.nii.gz"), tgdm.beforeReinitialization);
        System.out.println("Wrote isolated Java outer-surface output to " + output.getAbsolutePath());
    }

    private static final class CapturingTGDM extends MaGenericTGDM {
        float[][][] beforeReinitialization;

        CapturingTGDM() { super(ConnectivityRule.CONNECT_18_6); }

        @Override
        protected void endNarrowBand(ProfileTGDM profile) {
            float[][][] source = implicitLevelSet;
            beforeReinitialization = new float[source.length][source[0].length][source[0][0].length];
            for (int x = 0; x < source.length; x++)
                for (int y = 0; y < source[x].length; y++)
                    System.arraycopy(source[x][y], 0, beforeReinitialization[x][y], 0, source[x][y].length);
            super.endNarrowBand(profile);
        }
    }
}
