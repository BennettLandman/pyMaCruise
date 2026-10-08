import edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.KdPoint3;
import edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.KdTriangle;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

/** Compare selected JIST KdTriangle.distance(Point3f) on raw float32 meshes. */
public final class KdTriangleDistanceHarness {
    private static final class Triangle extends KdTriangle {
        Triangle(float[] v) {
            pts = new KdPoint3[] {
                new KdPoint3(v[0], v[1], v[2]),
                new KdPoint3(v[3], v[4], v[5]),
                new KdPoint3(v[6], v[7], v[8])
            };
            update();
        }

        @Override public void update() {
            minPoint.x = Math.min(Math.min(pts[0].x, pts[1].x), pts[2].x) - EPS;
            minPoint.y = Math.min(Math.min(pts[0].y, pts[1].y), pts[2].y) - EPS;
            minPoint.z = Math.min(Math.min(pts[0].z, pts[1].z), pts[2].z) - EPS;
            maxPoint.x = Math.max(Math.max(pts[0].x, pts[1].x), pts[2].x) + EPS;
            maxPoint.y = Math.max(Math.max(pts[0].y, pts[1].y), pts[2].y) + EPS;
            maxPoint.z = Math.max(Math.max(pts[0].z, pts[1].z), pts[2].z) + EPS;
        }
    }

    private static float[] readFloats(Path path) throws IOException {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length % 4 != 0) throw new IllegalArgumentException("float input is not aligned");
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[bytes.length / 4];
        for (int i = 0; i < values.length; i++) values[i] = buffer.getFloat();
        return values;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 4) throw new IllegalArgumentException(
            "usage: KdTriangleDistanceHarness triangles.raw points.raw face-count output.raw");
        float[] triangleValues = readFloats(Path.of(args[0]));
        float[] pointValues = readFloats(Path.of(args[1]));
        int faceCount = Integer.parseInt(args[2]);
        if (triangleValues.length != faceCount * 9 || pointValues.length % 3 != 0)
            throw new IllegalArgumentException("raw geometry count mismatch");
        ArrayList<Triangle> triangles = new ArrayList<>(faceCount);
        for (int i = 0; i < faceCount; i++) {
            float[] v = new float[9];
            System.arraycopy(triangleValues, i * 9, v, 0, 9);
            triangles.add(new Triangle(v));
        }
        ByteBuffer output = ByteBuffer.allocate((pointValues.length / 3) * 20)
            .order(ByteOrder.LITTLE_ENDIAN);
        for (int i = 0; i < pointValues.length; i += 3) {
            javax.vecmath.Point3f query = new javax.vecmath.Point3f(
                pointValues[i], pointValues[i + 1], pointValues[i + 2]);
            double bestDistance = Double.POSITIVE_INFINITY;
            javax.vecmath.Point3f bestPoint = null;
            for (Triangle triangle : triangles) {
                double distance = triangle.distance(query);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestPoint = KdTriangle.getLastIntersectionPoint();
                }
            }
            output.putDouble(bestDistance);
            output.putFloat(bestPoint.x);
            output.putFloat(bestPoint.y);
            output.putFloat(bestPoint.z);
        }
        Files.write(Path.of(args[3]), output.array());
    }
}
