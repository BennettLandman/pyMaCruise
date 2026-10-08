import edu.jhu.ece.iacl.algorithms.graphics.intersector.IntersectorTriangle;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import java.nio.*;
import java.nio.file.*;
import javax.vecmath.Point3f;

/** Compare selected Java KdTriangle distances and projected points by face id. */
public final class KdTriangleDistanceHarness {
    private static float[] readFloats(Path p) throws Exception {
        ByteBuffer b=ByteBuffer.wrap(Files.readAllBytes(p)).order(ByteOrder.LITTLE_ENDIAN);
        float[] a=new float[b.remaining()/4]; for(int i=0;i<a.length;i++) a[i]=b.getFloat(); return a;
    }
    private static int[] readInts(Path p) throws Exception {
        ByteBuffer b=ByteBuffer.wrap(Files.readAllBytes(p)).order(ByteOrder.LITTLE_ENDIAN);
        int[] a=new int[b.remaining()/4]; for(int i=0;i<a.length;i++) a[i]=b.getInt(); return a;
    }
    private static EmbeddedSurface mesh(Path v,Path f) throws Exception {
        float[] x=readFloats(v); int[] ids=readInts(f); Point3f[] p=new Point3f[x.length/3];
        for(int i=0;i<p.length;i++) p[i]=new Point3f(x[3*i],x[3*i+1],x[3*i+2]);
        return new EmbeddedSurface(p,ids);
    }
    public static void main(String[] a) throws Exception {
        EmbeddedSurface mesh=mesh(Path.of(a[0]),Path.of(a[1]));
        Point3f query=new Point3f(Float.parseFloat(a[2]),Float.parseFloat(a[3]),Float.parseFloat(a[4]));
        ByteBuffer out=ByteBuffer.allocate((a.length-6)*24).order(ByteOrder.LITTLE_ENDIAN);
        for(int i=6;i<a.length;i++) {
            int face=Integer.parseInt(a[i]); int at=3*face;
            IntersectorTriangle tri=new IntersectorTriangle(mesh.getCoordinateIndex(at),mesh.getCoordinateIndex(at+1),mesh.getCoordinateIndex(at+2),face,mesh);
            for(int warmup=0;warmup<12000;warmup++) tri.distance(query);
            double distance=tri.distance(query); Point3f p=edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.KdTriangle.getLastIntersectionPoint();
            out.putInt(face).putDouble(distance).putFloat(p.x).putFloat(p.y).putFloat(p.z);
        }
        Files.write(Path.of(a[5]),out.array());
    }
}
