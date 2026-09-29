/** Selected-method formula excerpt for LesionToadAtlas coefficient reduction. */
public final class LesionRegistrationVolumeAccumulationReferenceHarness {
    private static final int NX = 4, NY = 3, NZ = 5, CLASSES = 4, ND = 6;

    private static float membership(int x, int y, int z, int k) {
        float value = 0.21f;
        value += 0.09f * k;
        value += 0.017f * x;
        value += 0.011f * y;
        value += 0.007f * z;
        return value;
    }

    private static float prior(int x, int y, int z, int k) {
        if (k == 0) return 0.88f;
        if ((x + 2 * y + z + k) % 5 == 0) return 0.0f;
        float value = 0.11f;
        value += 0.04f * k;
        value += 0.003f * x;
        value -= 0.002f * y;
        value += 0.001f * z;
        return value;
    }

    private static float derivative(int axis, int x, int y, int z, int k) {
        if (axis == 0) return 0.001f * (k + x + 1);
        if (axis == 1) return 0.0013f * (k + y + 2);
        return 0.0007f * (k + z + 3);
    }

    private static float jacobian(int axis, int parameter, int x, int y, int z) {
        float value = 0.012f * (axis + 1) * (parameter + 1);
        value += 0.0001f * x;
        value -= 0.0002f * y;
        value += 0.0003f * z;
        return value;
    }

    private static String vector(float[] values) {
        StringBuilder out = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(Float.toString(values[i]));
        }
        return out.append(']').toString();
    }

    private static String run(String name, int offset, int step, boolean reverse) {
        float cost = 0.0f, norm = 0.0f;
        float[] gradient = new float[ND], hessian = new float[ND];
        int voxels = 0, terms = 0;
        int xStart = reverse ? NX - 1 - ((NX - 1 - offset) % step) : offset;
        int yStart = reverse ? NY - 1 - ((NY - 1 - offset) % step) : offset;
        int zStart = reverse ? NZ - 1 - ((NZ - 1 - offset) % step) : offset;
        int xDelta = reverse ? -step : step;
        int yDelta = reverse ? -step : step;
        int zDelta = reverse ? -step : step;
        for (int x = xStart; reverse ? x >= offset : x < NX; x += xDelta) {
            for (int y = yStart; reverse ? y >= offset : y < NY; y += yDelta) {
                for (int z = zStart; reverse ? z >= offset : z < NZ; z += zDelta) {
                    voxels++;
                    for (int k = reverse ? CLASSES - 1 : 1;
                         reverse ? k >= 1 : k < CLASSES;
                         k += reverse ? -1 : 1) {
                        float p = prior(x,y,z,k);
                        if (p > 0.0f) {
                            float[] dp = new float[ND];
                            for (int i=0; i<ND; i++) {
                                dp[i] = derivative(0,x,y,z,k)*jacobian(0,i,x,y,z)
                                    + derivative(1,x,y,z,k)*jacobian(1,i,x,y,z)
                                    + derivative(2,x,y,z,k)*jacobian(2,i,x,y,z);
                            }
                            float m = membership(x,y,z,k);
                            for (int i=0; i<ND; i++) {
                                hessian[i] += 2.0f*m*m*dp[i]*dp[i];
                                gradient[i] += 2.0f*m*m*p*dp[i];
                            }
                            cost += m*m*p*p;
                            norm += p*p;
                            terms++;
                        }
                    }
                }
            }
        }
        float[] normalizedGradient = new float[ND], normalizedHessian = new float[ND];
        for (int i=0; i<ND; i++) {
            normalizedGradient[i] = gradient[i] / norm;
            normalizedHessian[i] = hessian[i] / norm;
        }
        return "{\"name\":\""+name+"\",\"offset\":"+offset+",\"subsample\":"+step
            +",\"reverse\":"+reverse+",\"cost\":"+Float.toString(cost)
            +",\"norm\":"+Float.toString(norm)+",\"gradient\":"+vector(gradient)
            +",\"hessian\":"+vector(hessian)+",\"energy\":"+Float.toString(cost/norm)
            +",\"normalized_gradient\":"+vector(normalizedGradient)
            +",\"normalized_hessian\":"+vector(normalizedHessian)
            +",\"visited_voxels\":"+voxels+",\"contributing_class_terms\":"+terms+'}';
    }

    public static void main(String[] args) {
        System.out.println("{\"dimensions_xyz\":["+NX+","+NY+","+NZ+"],\"classes\":"+CLASSES
            +",\"cases\":["+run("dense",0,1,false)+","+run("subsampled",0,2,false)
            +","+run("reverse_control",0,1,true)+"]}");
    }
}
