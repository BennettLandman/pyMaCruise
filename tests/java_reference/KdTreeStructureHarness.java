package edu.jhu.ece.iacl.algorithms.graphics.intersector;

import edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.BBox;
import edu.jhu.ece.iacl.algorithms.graphics.locator.kdtree.KdTriangle;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import java.nio.*;
import java.nio.file.*;
import javax.vecmath.Point3f;

/** Dump the selected JIST tree's branch/leaf ordering for source-port checks. */
public final class KdTreeStructureHarness {
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
    private static void dump(BBox box,String path,StringBuilder out) {
        if(box instanceof IntersectorTriangle) {
            out.append(path).append("\tT\t").append(((IntersectorTriangle)box).fid).append('\n');
            return;
        }
        out.append(path).append("\tB\t").append(box.getDepth()).append('\t')
            .append(box.getChildren().size()).append('\t')
            .append(box.minPoint.x).append(',').append(box.minPoint.y).append(',').append(box.minPoint.z).append('\t')
            .append(box.maxPoint.x).append(',').append(box.maxPoint.y).append(',').append(box.maxPoint.z).append('\n');
        for(int i=0;i<box.getChildren().size();i++) dump(box.getChildren().get(i),path+i,out);
    }
    public static void main(String[] a) throws Exception {
        SurfaceIntersector si=new SurfaceIntersector(13,mesh(Path.of(a[0]),Path.of(a[1])));
        StringBuilder out=new StringBuilder();dump(si.tree.getRoot(),"r",out);
        Files.writeString(Path.of(a[2]),out.toString());
    }
}
