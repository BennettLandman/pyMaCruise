import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

/**
 * Dynamic selected-method excerpt for LesionToadPreprocess.uncropAndBuffer.
 *
 * The x/y/z copy loops and flattened index expression are copied from the
 * selected VandyCRUISE float, byte, and short overloads. Geometry and arrays
 * are supplied by the caller; this does not instantiate the MIPAV class or
 * serialize ModelImage/NIfTI files.
 */
public final class LesionOutputUncropBufferReferenceHarness {
    private static int nx, ny, nz, x0, xn, y0, yn, z0, zn, bs;

    private static int index(int x, int y, int z) {
        return x + nx * y + nx * ny * z;
    }

    private static void runFloat(DataInputStream in, DataOutputStream out, int channels) throws IOException {
        float[][][][] src = new float[channels][xn-x0+1+2*bs][yn-y0+1+2*bs][zn-z0+1+2*bs];
        for (int c=0;c<channels;c++)
            for (int x=0;x<src[c].length;x++)
                for (int y=0;y<src[c][x].length;y++)
                    for (int z=0;z<src[c][x][y].length;z++) src[c][x][y][z] = in.readFloat();
        float[] larger = new float[channels * nx * ny * nz];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
            for (int c=0;c<channels;c++) larger[index(x,y,z) + nx*ny*nz*c] = 0.0f;
        for (int x=Math.max(x0-bs,0);x<=Math.min(xn+bs,nx-1);x++)
            for (int y=Math.max(y0-bs,0);y<=Math.min(yn+bs,ny-1);y++)
                for (int z=Math.max(z0-bs,0);z<=Math.min(zn+bs,nz-1);z++)
                    for (int c=0;c<channels;c++)
                        larger[index(x,y,z) + nx*ny*nz*c] = src[c][x-x0+bs][y-y0+bs][z-z0+bs];
        for (float value : larger) out.writeFloat(value);
    }

    private static void runByte(DataInputStream in, DataOutputStream out) throws IOException {
        byte[][][] src = new byte[xn-x0+1+2*bs][yn-y0+1+2*bs][zn-z0+1+2*bs];
        for (int x=0;x<src.length;x++) for (int y=0;y<src[x].length;y++)
            for (int z=0;z<src[x][y].length;z++) src[x][y][z] = in.readByte();
        byte[] larger = new byte[nx*ny*nz];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) larger[index(x,y,z)] = 0;
        for (int x=Math.max(x0-bs,0);x<=Math.min(xn+bs,nx-1);x++)
            for (int y=Math.max(y0-bs,0);y<=Math.min(yn+bs,ny-1);y++)
                for (int z=Math.max(z0-bs,0);z<=Math.min(zn+bs,nz-1);z++)
                    larger[index(x,y,z)] = src[x-x0+bs][y-y0+bs][z-z0+bs];
        out.write(larger);
    }

    private static void runShort(DataInputStream in, DataOutputStream out) throws IOException {
        short[][][] src = new short[xn-x0+1+2*bs][yn-y0+1+2*bs][zn-z0+1+2*bs];
        for (int x=0;x<src.length;x++) for (int y=0;y<src[x].length;y++)
            for (int z=0;z<src[x][y].length;z++) src[x][y][z] = in.readShort();
        short[] larger = new short[nx*ny*nz];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) larger[index(x,y,z)] = 0;
        for (int x=Math.max(x0-bs,0);x<=Math.min(xn+bs,nx-1);x++)
            for (int y=Math.max(y0-bs,0);y<=Math.min(yn+bs,ny-1);y++)
                for (int z=Math.max(z0-bs,0);z<=Math.min(zn+bs,nz-1);z++)
                    larger[index(x,y,z)] = src[x-x0+bs][y-y0+bs][z-z0+bs];
        for (short value : larger) out.writeShort(value);
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 12) throw new IllegalArgumentException("type channels nx ny nz x0 xn y0 yn z0 zn padding");
        int type = Integer.parseInt(args[0]);
        int channels = Integer.parseInt(args[1]);
        nx = Integer.parseInt(args[2]); ny = Integer.parseInt(args[3]); nz = Integer.parseInt(args[4]);
        x0 = Integer.parseInt(args[5]); xn = Integer.parseInt(args[6]);
        y0 = Integer.parseInt(args[7]); yn = Integer.parseInt(args[8]);
        z0 = Integer.parseInt(args[9]); zn = Integer.parseInt(args[10]);
        bs = Integer.parseInt(args[11]);
        try (DataInputStream in = new DataInputStream(System.in); DataOutputStream out = new DataOutputStream(System.out)) {
            if (type == 0 && channels == 1) runByte(in, out);
            else if (type == 1 && channels == 1) runShort(in, out);
            else if (type == 2 && channels > 0) runFloat(in, out, channels);
            else throw new IllegalArgumentException("unsupported output type/channel count");
        }
    }
}
