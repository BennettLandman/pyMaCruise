/**
 * Pinned excerpts of LesionToadPreprocess.transformNewImage,
 * computeTemplateCoordinates, RotationMatrix.setParameters, and
 * ImageFunctions.linearInterpolation from the selected VandyCRUISE source.
 */
public final class LesionToadsTransformNewImageReferenceHarness {
    private static float[][] rotation(float a, float b, float c) {
        double norm = Math.sqrt(a * a + b * b + c * c);
        double an = a, bn = b, cn = c;
        if (norm > 1.0) {
            double m = norm % 2;
            if (m > 1.0) {
                an = m * a / norm - 2 * a / norm;
                bn = m * b / norm - 2 * b / norm;
                cn = m * c / norm - 2 * c / norm;
            } else {
                an = a * m / norm; bn = b * m / norm; cn = c * m / norm;
            }
        }
        double a2 = an * an, b2 = bn * bn, c2 = cn * cn;
        double w2 = 1.0 - a2 - b2 - c2;
        if (w2 < 0.0 && w2 > -1e-10) w2 = 0.0;
        double w = Math.sqrt(w2);
        return new float[][] {
            {(float)(1 - 2*b2 - 2*c2), (float)(2*an*bn - 2*w*cn), (float)(2*an*cn + 2*w*bn)},
            {(float)(2*an*bn + 2*w*cn), (float)(1 - 2*a2 - 2*c2), (float)(2*bn*cn - 2*w*an)},
            {(float)(2*an*cn - 2*w*bn), (float)(2*bn*cn + 2*w*an), (float)(1 - 2*a2 - 2*b2)}
        };
    }

    private static float linear(float[][][] image, float value, float x, float y, float z) {
        int nx = image.length, ny = image[0].length, nz = image[0][0].length;
        if (x < 0 || x > nx-2 || y < 0 || y > ny-2 || z < 0 || z > nz-2) return value;
        int x0 = (int)Math.floor(x), y0 = (int)Math.floor(y), z0 = (int)Math.floor(z);
        float alpha = x - x0, nalpha = 1.0f - alpha;
        float beta = y - y0, nbeta = 1.0f - beta;
        float gamma = z - z0, ngamma = 1.0f - gamma;
        float val = nalpha*nbeta*ngamma*image[x0][y0][z0]
            + alpha*nbeta*ngamma*image[x0+1][y0][z0]
            + nalpha*beta*ngamma*image[x0][y0+1][z0]
            + nalpha*nbeta*gamma*image[x0][y0][z0+1]
            + alpha*beta*ngamma*image[x0+1][y0+1][z0]
            + nalpha*beta*gamma*image[x0][y0+1][z0+1]
            + alpha*nbeta*gamma*image[x0+1][y0][z0+1]
            + alpha*beta*gamma*image[x0+1][y0+1][z0+1];
        return val;
    }

    private static float[][][] transformNewImage(float[][][] image, float background,
            float orx, float ory, float orz, float[] transform) {
        int nix = 7, niy = 8, niz = 9;
        int ntx = 9, nty = 10, ntz = 11;
        int x0i = 1, xNi = 5, y0i = 1, yNi = 6, z0i = 1, zNi = 7, bs = 1;
        float rix = 1.1f, riy = 1.3f, riz = 0.9f;
        float rtx = 0.9f, rty = 1.0f, rtz = 1.1f;
        float xCi = nix/2.0f, yCi = niy/2.0f, zCi = niz/2.0f;
        float xCt = ntx/2.0f, yCt = nty/2.0f, zCt = ntz/2.0f;
        float[][] rot = rotation(transform[0], transform[1], transform[2]);
        int mx=xNi-x0i+1+2*bs, my=yNi-y0i+1+2*bs, mz=zNi-z0i+1+2*bs;
        float[][][] smaller = new float[mx][my][mz];
        for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) smaller[x][y][z]=background;
        for (int x=x0i;x<=xNi;x++) for (int y=y0i;y<=yNi;y++) for (int z=z0i;z<=zNi;z++) {
            float xt = (rot[0][0]*(x-xCi)*rix/orx + rot[0][1]*(y-yCi)*riy/ory + rot[0][2]*(z-zCi)*riz/orz + transform[3])/rtx + xCt;
            float yt = (rot[1][0]*(x-xCi)*rix/orx + rot[1][1]*(y-yCi)*riy/ory + rot[1][2]*(z-zCi)*riz/orz + transform[4])/rty + yCt;
            float zt = (rot[2][0]*(x-xCi)*rix/orx + rot[2][1]*(y-yCi)*riy/ory + rot[2][2]*(z-zCi)*riz/orz + transform[5])/rtz + zCt;
            smaller[x-x0i+bs][y-y0i+bs][z-z0i+bs] = linear(image, background, xt, yt, zt);
        }
        return smaller;
    }

    public static void main(String[] args) {
        float[][][] image = new float[9][10][11];
        for (int x=0;x<9;x++) for (int y=0;y<10;y++) for (int z=0;z<11;z++)
            image[x][y][z] = ((x*7 + y*3 + z*5 + x*y + z*y*2) % 97) / 16.0f;
        float[][][] output = transformNewImage(image, -2.5f, 1.0f, 1.0f, 1.0f,
            new float[]{0.02f, -0.03f, 0.01f, 0.3f, -0.15f, 0.4f});
        float[][][] negativeX = transformNewImage(image, -2.5f, -1.0f, 1.0f, 1.0f,
            new float[]{0.02f, -0.03f, 0.01f, 0.3f, -0.15f, 0.4f});
        float[][][] outOfDomain = transformNewImage(image, -2.5f, 1.0f, 1.0f, 1.0f,
            new float[]{0.02f, -0.03f, 0.01f, 2.0f, -0.15f, 0.4f});
        System.out.print("{\"source_method\":\"LesionToadPreprocess.transformNewImage\",");
        System.out.print("\"source_sha256\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",");
        System.out.print("\"image_functions_sha256\":\"be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892\",");
        System.out.print("\"rotation_matrix_sha256\":\"144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0\",");
        System.out.print("\"template_shape\":[9,10,11],\"image_shape\":[7,8,9],\"bounds\":[1,5,1,6,1,7],\"padding\":1,\"background\":-2.5,\"output_shape\":[7,8,9],\"output_bits\":[");
        boolean first = true;
        for (int x=0;x<output.length;x++) for (int y=0;y<output[0].length;y++) for (int z=0;z<output[0][0].length;z++) {
            if (!first) System.out.print(","); first=false;
            System.out.print(Integer.toUnsignedLong(Float.floatToIntBits(output[x][y][z])));
        }
        System.out.print("],\"negative_x_orientation_bits\":[");
        first = true;
        for (int x=0;x<negativeX.length;x++) for (int y=0;y<negativeX[0].length;y++) for (int z=0;z<negativeX[0][0].length;z++) {
            if (!first) System.out.print(","); first=false;
            System.out.print(Integer.toUnsignedLong(Float.floatToIntBits(negativeX[x][y][z])));
        }
        System.out.print("],\"out_of_domain_bits\":[");
        first = true;
        for (int x=0;x<outOfDomain.length;x++) for (int y=0;y<outOfDomain[0].length;y++) for (int z=0;z<outOfDomain[0][0].length;z++) {
            if (!first) System.out.print(","); first=false;
            System.out.print(Integer.toUnsignedLong(Float.floatToIntBits(outOfDomain[x][y][z])));
        }
        System.out.println("]}");
    }
}
