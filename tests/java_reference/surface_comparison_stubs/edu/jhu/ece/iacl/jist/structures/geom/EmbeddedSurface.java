package edu.jhu.ece.iacl.jist.structures.geom;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

/** Capture the selected extractor's returned vertex and face arrays verbatim. */
public final class EmbeddedSurface {
    public enum Direction { COUNTER_CLOCKWISE, CLOCKWISE }

    public final Point3f[] points;
    public final int[] indices;

    public EmbeddedSurface(Point3f[] points, Vector3f[] normals,
                           Color3f[] colors, int[] indices) {
        this.points = points;
        this.indices = indices;
    }

    public void setName(String name) {}
}
