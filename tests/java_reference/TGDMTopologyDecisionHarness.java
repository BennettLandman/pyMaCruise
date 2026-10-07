import java.io.File;

import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.vanderbilt.masi.plugins.CRUISE.tgdm.MaGenericTGDM;

/** Expose the selected TGDM simple-point decision for one supplied sign field. */
public final class TGDMTopologyDecisionHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("usage: TGDMTopologyDecisionHarness CENTRAL_NIFTI");
        float[][][] values = TGDMCoreReferenceHarness.readNifti(new File(args[0]));
        int x = 73, y = 76, z = 73;
        byte[][][] status = new byte[values.length][values[0].length][values[0][0].length];
        for (int i = 0; i < values.length; i++)
            for (int j = 0; j < values[i].length; j++)
                for (int k = 0; k < values[i][j].length; k++)
                    status[i][j][k] = (byte) (values[i][j][k] > 0 ? -1 : 1);
        ExposedTGDM tgdm = new ExposedTGDM();
        System.out.println("java_simple_point=" + tgdm.isSimple(x, y, z, status));
    }

    private static final class ExposedTGDM extends MaGenericTGDM {
        ExposedTGDM() { super(ConnectivityRule.CONNECT_18_6); }

        boolean isSimple(int x, int y, int z, byte[][][] status) {
            int[][][] cube = new int[3][3][3];
            smpTest.test(x, y, z, cube, status);
            int foreground = fgRule.check(cube, 1);
            int background = bgRule.check(cube, -1);
            System.out.println("foreground_components=" + foreground);
            System.out.println("background_components=" + background);
            return foreground == 1 && background == 1;
        }
    }
}
