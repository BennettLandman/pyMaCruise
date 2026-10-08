import edu.jhu.ece.iacl.algorithms.graphics.intersector.IntersectorTriangle;
import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import java.nio.*;
import java.nio.file.*;
import javax.vecmath.Point3f;

/** Return selected JIST triangle ids for a batch of closest-point queries. */
public final class SurfaceIntersectorFaceQueryHarness {
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
        EmbeddedSurface query=mesh(Path.of(a[0]),Path.of(a[1]));
        EmbeddedSurface target=mesh(Path.of(a[2]),Path.of(a[3]));
        SurfaceIntersector si=new SurfaceIntersector(13,target);
        ByteBuffer out=ByteBuffer.allocate(query.getVertexCount()*8).order(ByteOrder.LITTLE_ENDIAN);
        for(int i=0;i<query.getVertexCount();i++) {
            double distance=si.distance(query.getVertex(i));
            IntersectorTriangle triangle=(IntersectorTriangle)si.getLastIntersectionTriangle();
            out.putInt(triangle.fid).putFloat((float)distance);
        }
        Files.write(Path.of(a[4]),out.array());
    }
}
