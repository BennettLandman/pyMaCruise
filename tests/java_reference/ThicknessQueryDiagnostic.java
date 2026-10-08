import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.algorithms.graphics.intersector.IntersectorTriangle;
import edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.KdTriangle;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import java.nio.*;
import java.nio.file.*;
import javax.vecmath.Point3f;

/** Expose selected-source closest-point intermediates for one thickness query. */
public final class ThicknessQueryDiagnostic {
    private static float[] floats(Path p) throws Exception {
        ByteBuffer b=ByteBuffer.wrap(Files.readAllBytes(p)).order(ByteOrder.LITTLE_ENDIAN);
        float[] a=new float[b.remaining()/4]; for(int i=0;i<a.length;i++) a[i]=b.getFloat(); return a;
    }
    private static int[] ints(Path p) throws Exception {
        ByteBuffer b=ByteBuffer.wrap(Files.readAllBytes(p)).order(ByteOrder.LITTLE_ENDIAN);
        int[] a=new int[b.remaining()/4]; for(int i=0;i<a.length;i++) a[i]=b.getInt(); return a;
    }
    private static EmbeddedSurface mesh(Path v,Path f) throws Exception {
        float[] x=floats(v); int[] ids=ints(f); Point3f[] p=new Point3f[x.length/3];
        for(int i=0;i<p.length;i++) p[i]=new Point3f(x[3*i],x[3*i+1],x[3*i+2]);
        return new EmbeddedSurface(p,ids);
    }
    private static void hit(ByteBuffer b,SurfaceIntersector si,double d) {
        b.putDouble(d); Point3f p=si.getLastIntersectionPoint();
        b.putFloat(p.x).putFloat(p.y).putFloat(p.z);
        IntersectorTriangle t=(IntersectorTriangle)si.getLastIntersectionTriangle();
        b.putInt(t.fid).putInt(t.id1).putInt(t.id2).putInt(t.id3);
    }
    public static void main(String[] a) throws Exception {
        EmbeddedSurface inner=mesh(Path.of(a[0]),Path.of(a[1]));
        EmbeddedSurface outer=mesh(Path.of(a[2]),Path.of(a[3])); int index=Integer.parseInt(a[4]);
        Point3f query=inner.getVertex(index);
        SurfaceIntersector toOuter=new SurfaceIntersector((AbstractCalculation)null,13,outer);
        SurfaceIntersector toInner=new SurfaceIntersector((AbstractCalculation)null,13,inner);
        double d1=toOuter.distance(query); Point3f outerHit=toOuter.getLastIntersectionPoint();
        double d2=toInner.distance(outerHit);
        ByteBuffer b=ByteBuffer.allocate(12+2*(8+12+16)).order(ByteOrder.LITTLE_ENDIAN);
        b.putFloat(query.x).putFloat(query.y).putFloat(query.z); hit(b,toOuter,d1); hit(b,toInner,d2);
        Files.write(Path.of(a[5]),b.array());
    }
}
