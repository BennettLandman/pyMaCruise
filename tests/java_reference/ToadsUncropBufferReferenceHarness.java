/*
 * Selected-method excerpt harness for VandyCRUISE MultipleToadPreprocess.
 * The three uncropAndBuffer overload bodies below are copied from the
 * selected source file (lines 478-535), with fields passed as arguments.
 * This intentionally excludes MIPAV/JIST serialization and orchestration.
 */
public final class ToadsUncropBufferReferenceHarness {
    static final int NX = 4, NY = 5, NZ = 3;
    static final int X0 = -1, XN = 1, Y0 = 1, YN = 2, Z0 = -1, ZN = 0, BS = 1;

    static float[] uncrop(float[][][] src, float bgVal) {
        float[] larger = new float[NX * NY * NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            larger[x + NX*y + NX*NY*z] = bgVal;
        for (int x=Math.max(X0-BS,0);x<=Math.min(XN+BS,NX-1);x++)
            for (int y=Math.max(Y0-BS,0);y<=Math.min(YN+BS,NY-1);y++)
                for (int z=Math.max(Z0-BS,0);z<=Math.min(ZN+BS,NZ-1);z++)
                    larger[x + NX*y + NX*NY*z] = src[x-X0+BS][y-Y0+BS][z-Z0+BS];
        return larger;
    }

    static byte[] uncrop(byte[][][] src, byte bgVal) {
        byte[] larger = new byte[NX * NY * NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            larger[x + NX*y + NX*NY*z] = bgVal;
        for (int x=Math.max(X0-BS,0);x<=Math.min(XN+BS,NX-1);x++)
            for (int y=Math.max(Y0-BS,0);y<=Math.min(YN+BS,NY-1);y++)
                for (int z=Math.max(Z0-BS,0);z<=Math.min(ZN+BS,NZ-1);z++)
                    larger[x + NX*y + NX*NY*z] = src[x-X0+BS][y-Y0+BS][z-Z0+BS];
        return larger;
    }

    static float[] uncrop(float[][][][] src, float bgVal, int nv) {
        float[] larger = new float[nv * NX * NY * NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) for (int v=0;v<nv;v++)
            larger[x + NX*y + NX*NY*z + NX*NY*NZ*v] = bgVal;
        for (int x=Math.max(X0-BS,0);x<=Math.min(XN+BS,NX-1);x++)
            for (int y=Math.max(Y0-BS,0);y<=Math.min(YN+BS,NY-1);y++)
                for (int z=Math.max(Z0-BS,0);z<=Math.min(ZN+BS,NZ-1);z++)
                    for (int v=0;v<nv;v++)
                        larger[x + NX*y + NX*NY*z + NX*NY*NZ*v] = src[v][x-X0+BS][y-Y0+BS][z-Z0+BS];
        return larger;
    }

    static void floats(String name, float[] values) {
        System.out.print("\"" + name + "\":[");
        for (int i=0;i<values.length;i++) { if (i>0) System.out.print(","); System.out.print(Float.toString(values[i])); }
        System.out.print("]");
    }
    static void bytes(String name, byte[] values) {
        System.out.print("\"" + name + "\":[");
        for (int i=0;i<values.length;i++) { if (i>0) System.out.print(","); System.out.print(values[i]); }
        System.out.print("]");
    }

    public static void main(String[] args) {
        int cx=XN-X0+1+2*BS, cy=YN-Y0+1+2*BS, cz=ZN-Z0+1+2*BS;
        float[][][] scalar = new float[cx][cy][cz];
        byte[][][] labels = new byte[cx][cy][cz];
        float[][][][] channels = new float[2][cx][cy][cz];
        for (int x=0;x<cx;x++) for (int y=0;y<cy;y++) for (int z=0;z<cz;z++) {
            scalar[x][y][z] = x*100.0f + y*10.0f + z + 0.25f;
            labels[x][y][z] = (byte)(100 + 20*x + 5*y + z);
            channels[0][x][y][z] = x*1000.0f + y*100.0f + z + 0.5f;
            channels[1][x][y][z] = -x*1000.0f - y*100.0f - z - 0.75f;
        }
        System.out.print("{\"dims\":["+NX+","+NY+","+NZ+"],\"bounds\":["+X0+","+XN+","+Y0+","+YN+","+Z0+","+ZN+"],\"padding\":"+BS+",\"crop_dims\":["+cx+","+cy+","+cz+"],");
        floats("scalar", uncrop(scalar, -7.5f)); System.out.print(",");
        bytes("labels", uncrop(labels, (byte)9)); System.out.print(",");
        floats("channels", uncrop(channels, 2.25f, 2));
        System.out.println("}");
    }
}
