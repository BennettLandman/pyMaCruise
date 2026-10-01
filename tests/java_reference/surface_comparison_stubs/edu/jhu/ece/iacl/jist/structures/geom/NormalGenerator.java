package edu.jhu.ece.iacl.jist.structures.geom;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

/** Normals are not used in the geometry comparison; preserve the constructor contract. */
public final class NormalGenerator {
    public static Vector3f[] generate(Point3f[] points, int[] indices) {
        Vector3f[] normals = new Vector3f[points.length];
        for (int i = 0; i < normals.length; i++) normals[i] = new Vector3f();
        return normals;
    }
}
