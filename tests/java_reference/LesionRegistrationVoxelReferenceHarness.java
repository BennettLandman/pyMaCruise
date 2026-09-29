/*
 * Test-only source-formula excerpt for one voxel in the selected
 * LesionToadAtlas.computeRegistrationCoefficients loop. It reproduces the
 * class-zero exclusion, positive-prior gate, prior-square cost/norm, and
 * diagonal gradient/Hessian updates for identity rotation at the image center.
 */
public final class LesionRegistrationVoxelReferenceHarness {
    private static final int N = 6;

    private static float value(int clazz, int x, int y, int z) {
        if (clazz == 0) return 0.9f;
        if (clazz == 1) {
            return 0.125f + 0.05f*x + 0.025f*y + 0.0125f*z;
        }
        if (clazz == 2) {
            return 0.1f + 0.025f*x + 0.05f*y + 0.0125f*z;
        }
        return 0.0f;
    }

    private static float dx(int clazz, int x, int y, int z) {
        if (x > N-2 || y > N-2 || z > N-2) return 0.0f;
        return value(clazz, x+1, y, z) - value(clazz, x, y, z);
    }

    private static float dy(int clazz, int x, int y, int z) {
        if (x > N-2 || y > N-2 || z > N-2) return 0.0f;
        return value(clazz, x, y+1, z) - value(clazz, x, y, z);
    }

    private static float dz(int clazz, int x, int y, int z) {
        if (x > N-2 || y > N-2 || z > N-2) return 0.0f;
        return value(clazz, x, y, z+1) - value(clazz, x, y, z);
    }

    public static void main(String[] args) {
        final int x=3, y=3, z=3;
        float[] memberships={0.9f,0.5f,0.75f,0.6f};
        float cost=0.0f, norm=0.0f;
        float[] gradient=new float[6], hessian=new float[6];
        for (int k=1; k<memberships.length; k++) {
            float prior=value(k,x,y,z);
            if (prior>0.0f) {
                float[] dprior={dx(k,x,y,z),dy(k,x,y,z),dz(k,x,y,z)};
                float[] dp={0,0,0,dprior[0],dprior[1],dprior[2]};
                float m=memberships[k];
                cost += m*m*prior*prior;
                norm += prior*prior;
                for (int i=0; i<6; i++) {
                    gradient[i] += 2.0f*m*m*prior*dp[i];
                    hessian[i] += 2.0f*m*m*dp[i]*dp[i];
                }
            }
        }
        StringBuilder out=new StringBuilder("{\"priors\":[");
        for (int k=0; k<memberships.length; k++) {
            if (k>0) out.append(',');
            out.append(Float.toString(k==0?0.0f:value(k,x,y,z)));
        }
        out.append("],\"prior_gradients_xyz\":[[0,0,0]");
        for (int k=1; k<memberships.length; k++) {
            out.append(",[").append(Float.toString(dx(k,x,y,z))).append(',')
               .append(Float.toString(dy(k,x,y,z))).append(',')
               .append(Float.toString(dz(k,x,y,z))).append(']');
        }
        out.append("],\"cost\":").append(Float.toString(cost))
            .append(",\"norm\":").append(Float.toString(norm));
        for (String name : new String[]{"gradient","hessian"}) {
            float[] values=name.equals("gradient")?gradient:hessian;
            out.append(",\"").append(name).append("\":[");
            for (int i=0; i<values.length; i++) {
                if (i>0) out.append(',');
                out.append(Float.toString(values[i]));
            }
            out.append(']');
        }
        System.out.println(out.append('}'));
    }
}
