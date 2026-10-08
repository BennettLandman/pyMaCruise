import edu.vanderbilt.masi.algorithms.CRUISE.utilities.ClosestPointThicknessForInner;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import javax.vecmath.Point3f;

/** Execute the selected VandyCRUISE ClosestPointThicknessForInner mesh overload. */
public final class ClosestPointThicknessForInnerHarness {
    private static float[] readFloats(Path path) throws Exception {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length % 4 != 0) throw new IllegalArgumentException("float input is not aligned");
        ByteBuffer in = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[bytes.length / 4];
        for (int i = 0; i < values.length; i++) values[i] = in.getFloat();
        return values;
    }
    private static int[] readInts(Path path) throws Exception {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length % 4 != 0) throw new IllegalArgumentException("index input is not aligned");
        ByteBuffer in = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        int[] values = new int[bytes.length / 4];
        for (int i = 0; i < values.length; i++) values[i] = in.getInt();
        return values;
    }
    private static EmbeddedSurface mesh(Path vertexPath, Path facePath) throws Exception {
        float[] raw = readFloats(vertexPath);
        int[] faces = readInts(facePath);
        if (raw.length % 3 != 0 || faces.length % 3 != 0)
            throw new IllegalArgumentException("mesh geometry count mismatch");
        Point3f[] points = new Point3f[raw.length / 3];
        for (int i = 0; i < points.length; i++)
            points[i] = new Point3f(raw[3*i], raw[3*i+1], raw[3*i+2]);
        return new EmbeddedSurface(points, faces);
    }
    public static void main(String[] args) throws Exception {
        if (args.length != 7) throw new IllegalArgumentException(
            "usage: harness inner-v inner-f central-v central-f outer-v outer-f out.raw");
        EmbeddedSurface inner = mesh(Path.of(args[0]), Path.of(args[1]));
        EmbeddedSurface central = mesh(Path.of(args[2]), Path.of(args[3]));
        EmbeddedSurface outer = mesh(Path.of(args[4]), Path.of(args[5]));
        ClosestPointThicknessForInner solver = new ClosestPointThicknessForInner();
        solver.solve(inner, central, outer);
        inner = solver.getInnerSurface();
        central = solver.getCentralSurface();
        int innerCount = inner.getVertexCount();
        int centralCount = central.getVertexCount();
        ByteBuffer out = ByteBuffer.allocate(8 + innerCount * 24 + centralCount * 32)
            .order(ByteOrder.LITTLE_ENDIAN);
        out.putInt(innerCount);
        for (int i = 0; i < innerCount; i++)
            for (double value : inner.getVertexData(i)) out.putDouble(value);
        out.putInt(centralCount);
        for (int i = 0; i < centralCount; i++)
            for (double value : central.getVertexData(i)) out.putDouble(value);
        Files.write(Path.of(args[6]), out.array());
    }
}
