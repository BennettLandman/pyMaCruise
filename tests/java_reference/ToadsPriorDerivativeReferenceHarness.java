/*
 * Test-only excerpt of ImageFunctions.linearInterpolationXderivative,
 * linearInterpolationYderivative, and linearInterpolationZderivative for the
 * selected VandyCRUISE atlas registration. The formulas below retain the
 * source branch order and float32 evaluation, including zero derivatives
 * outside [0, size-2] on any axis.
 */
public final class ToadsPriorDerivativeReferenceHarness {
    private static final int N = 5;
    private static final float[][][] shape = new float[N][N][N];

    private static float dx(float x, float y, float z) {
        if (x < 0 || x > N-2 || y < 0 || y > N-2 || z < 0 || z > N-2) return 0.0f;
        int x0=(int)Math.floor(x), y0=(int)Math.floor(y), z0=(int)Math.floor(z);
        float beta=y-y0, gamma=z-z0, nbeta=1.0f-beta, ngamma=1.0f-gamma;
        return nbeta*ngamma*(shape[x0+1][y0][z0]-shape[x0][y0][z0])
             + beta*ngamma*(shape[x0+1][y0+1][z0]-shape[x0][y0+1][z0])
             + nbeta*gamma*(shape[x0+1][y0][z0+1]-shape[x0][y0][z0+1])
             + beta*gamma*(shape[x0+1][y0+1][z0+1]-shape[x0][y0+1][z0+1]);
    }

    private static float dy(float x, float y, float z) {
        if (x < 0 || x > N-2 || y < 0 || y > N-2 || z < 0 || z > N-2) return 0.0f;
        int x0=(int)Math.floor(x), y0=(int)Math.floor(y), z0=(int)Math.floor(z);
        float alpha=x-x0, gamma=z-z0, nalpha=1.0f-alpha, ngamma=1.0f-gamma;
        return nalpha*ngamma*(shape[x0][y0+1][z0]-shape[x0][y0][z0])
             + alpha*ngamma*(shape[x0+1][y0+1][z0]-shape[x0+1][y0][z0])
             + nalpha*gamma*(shape[x0][y0+1][z0+1]-shape[x0][y0][z0+1])
             + alpha*gamma*(shape[x0+1][y0+1][z0+1]-shape[x0+1][y0][z0+1]);
    }

    private static float dz(float x, float y, float z) {
        if (x < 0 || x > N-2 || y < 0 || y > N-2 || z < 0 || z > N-2) return 0.0f;
        int x0=(int)Math.floor(x), y0=(int)Math.floor(y), z0=(int)Math.floor(z);
        float alpha=x-x0, beta=y-y0, nalpha=1.0f-alpha, nbeta=1.0f-beta;
        return nalpha*nbeta*(shape[x0][y0][z0+1]-shape[x0][y0][z0])
             + alpha*nbeta*(shape[x0+1][y0][z0+1]-shape[x0+1][y0][z0])
             + nalpha*beta*(shape[x0][y0+1][z0+1]-shape[x0][y0+1][z0])
             + alpha*beta*(shape[x0+1][y0+1][z0+1]-shape[x0+1][y0+1][z0]);
    }

    public static void main(String[] args) {
        for (int x=0; x<N; x++) for (int y=0; y<N; y++) for (int z=0; z<N; z++)
            shape[x][y][z] = ((7*x + 3*y + 5*z + 2*x*y + 4*y*z) % 23) / 22.0f;
        float[][] points = {{0.25f,1.5f,2.75f}, {3.0f,3.0f,3.0f}, {3.25f,1.0f,1.0f}, {-0.1f,1.0f,1.0f}};
        StringBuilder out = new StringBuilder("{\"derivatives\":[");
        for (int i=0; i<points.length; i++) {
            if (i>0) out.append(',');
            float[] p=points[i];
            out.append('[').append(Float.toString(dx(p[0],p[1],p[2]))).append(',')
               .append(Float.toString(dy(p[0],p[1],p[2]))).append(',')
               .append(Float.toString(dz(p[0],p[1],p[2]))).append(']');
        }
        System.out.println(out.append("]}"));
    }
}
