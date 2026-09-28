import edu.jhu.ece.iacl.algorithms.tgdm.DynamicPressureCortexProfile;
import edu.jhu.ece.iacl.algorithms.tgdm.DynamicPressureTGDM;
import edu.jhu.ece.iacl.algorithms.tgdm.GenericTGDM.NBPoint;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/**
 * Runs the selected VandyCRUISE DynamicPressureTGDM.updatePoint source on
 * deterministic scalar neighborhoods. The superclass seam supplies state;
 * this harness does not claim to run GenericTGDM's evolution loop.
 */
public final class DynamicPressureUpdateReferenceHarness {
    private static final int N = 5;

    private static final class Probe extends DynamicPressureTGDM {
        Probe() { super(ConnectivityRule.CONNECT_18_6); }

        double update(float[][][] phi, float[][][] pressure,
                float[][][] reference, int x, int y, int z,
                double dt, double alpha, double alignment, double curvature) {
            rows = phi.length;
            cols = phi[0].length;
            slices = phi[0][0].length;
            implicitLevelSet = phi;
            regionForces = pressure;
            timeStep = dt;
            this.alpha = (float) alpha;
            alignmentForce = (float) alignment;
            setReferenceLevelSet(new ImageDataFloat(reference));
            NBPoint point = new NBPoint();
            point.x = x;
            point.y = y;
            point.z = z;
            return updatePoint(point,
                    new DynamicPressureCortexProfile(-curvature, 0.0, 0.5, 5));
        }
    }

    public static void main(String[] args) {
        float[][][] phi = new float[N][N][N];
        float[][][] reference = new float[N][N][N];
        float[][][] pressurePositive = new float[N][N][N];
        float[][][] pressureNegative = new float[N][N][N];
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    double dx = x - 2.0, dy = y - 2.0, dz = z - 2.0;
                    phi[x][y][z] = (float) (0.18 * dx + 0.07 * dy * dy + 0.03 * dy * dz);
                    reference[x][y][z] = (float) (0.15 * dx + 0.12 * dz + 0.03 * dy * dy);
                    pressurePositive[x][y][z] = (float) (0.35 + 0.12 * dx - 0.08 * dy + 0.04 * dz);
                    pressureNegative[x][y][z] = (float) (-0.35 + 0.12 * dx - 0.08 * dy + 0.04 * dz);
                }
            }
        }
        Probe probe = new Probe();
        double center = probe.update(phi, pressurePositive, reference,
                2, 2, 2, 0.1, 0.5, 1.0, 0.5);
        double boundary = probe.update(phi, pressureNegative, reference,
                0, 2, 2, 0.1, 0.5, 1.0, 0.5);
        System.out.printf(
                "{\"source\":\"VandyCRUISE DynamicPressureTGDM.updatePoint\","
                + "\"source_revision\":\"1.1\",\"shape\":[5,5,5],"
                + "\"phantom\":{\"phi\":\"0.18*dx+0.07*dy^2+0.03*dy*dz\","
                + "\"reference\":\"0.15*dx+0.12*dz+0.03*dy^2\","
                + "\"pressure\":\"offset+0.12*dx-0.08*dy+0.04*dz\","
                + "\"center_offset\":0.35,\"boundary_offset\":-0.35},"
                + "\"parameters\":{\"time_step\":0.1,\"alpha\":0.5,"
                + "\"alignment_force\":1.0,\"curvature_profile_force\":-0.5,"
                + "\"pressure_profile_force\":0.5},"
                + "\"results\":[{\"name\":\"interior_positive_pressure\","
                + "\"xyz\":[2,2,2],\"old_phi\":%.9g,\"region_force\":%.9g,\"updated_phi\":%.17g},"
                + "{\"name\":\"boundary_negative_pressure\","
                + "\"xyz\":[0,2,2],\"old_phi\":%.9g,\"region_force\":%.9g,\"updated_phi\":%.17g}]}",
                phi[2][2][2], pressurePositive[2][2][2], center,
                phi[0][2][2], pressureNegative[0][2][2], boundary);
    }
}
