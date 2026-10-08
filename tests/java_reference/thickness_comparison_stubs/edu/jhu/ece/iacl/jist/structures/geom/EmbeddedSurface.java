package edu.jhu.ece.iacl.jist.structures.geom;

import javax.vecmath.Point3f;

/** Minimal coordinate/index shell for selected SurfaceIntersector source. */
public final class EmbeddedSurface {
    private final Point3f[] points;
    private final int[] indices;

    public EmbeddedSurface(Point3f[] points, int[] indices) {
        this.points = points;
        this.indices = indices;
    }

    public int getIndexCount() { return indices.length; }
    public int getCoordinateIndex(int index) { return indices[index]; }
    public Point3f getVertex(int index) { return points[index]; }
}
