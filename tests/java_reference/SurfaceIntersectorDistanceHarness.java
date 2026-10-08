import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import javax.vecmath.Point3f;

/** Invoke selected SurfaceIntersector.distance on raw float32 geometry. */
public final class SurfaceIntersectorDistanceHarness {
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

    public static void main(String[] args) throws Exception {
        if (args.length != 4) throw new IllegalArgumentException(
            "usage: SurfaceIntersectorDistanceHarness vertices.raw faces.raw queries.raw output.raw");
        float[] rawVertices = readFloats(Path.of(args[0]));
        int[] faces = readInts(Path.of(args[1]));
        float[] rawQueries = readFloats(Path.of(args[2]));
        if (rawVertices.length % 3 != 0 || faces.length % 3 != 0 || rawQueries.length % 3 != 0)
            throw new IllegalArgumentException("geometry count mismatch");
        Point3f[] vertices = new Point3f[rawVertices.length / 3];
        for (int i = 0; i < vertices.length; i++)
            vertices[i] = new Point3f(rawVertices[3*i], rawVertices[3*i+1], rawVertices[3*i+2]);
        EmbeddedSurface mesh = new EmbeddedSurface(vertices, faces);
        SurfaceIntersector intersector = new SurfaceIntersector(13, mesh);
        ByteBuffer output = ByteBuffer.allocate((rawQueries.length / 3) * 20)
            .order(ByteOrder.LITTLE_ENDIAN);
        for (int i = 0; i < rawQueries.length; i += 3) {
            Point3f query = new Point3f(rawQueries[i], rawQueries[i+1], rawQueries[i+2]);
            double distance = intersector.distance(query);
            Point3f point = intersector.getLastIntersectionPoint();
            output.putDouble(distance);
            output.putFloat(point.x);
            output.putFloat(point.y);
            output.putFloat(point.z);
        }
        Files.write(Path.of(args[3]), output.array());
    }
}
