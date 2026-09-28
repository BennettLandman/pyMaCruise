package edu.jhu.ece.iacl.algorithms.graphics;

import javax.vecmath.Vector3d;
import javax.vecmath.Point3f;

/** Numerical helper subset used by the unmodified VandyCRUISE profiles. */
public final class GeometricUtilities {
    private GeometricUtilities() {}
    public static void normalize(Vector3d value) {
        double length = value.length();
        if (length > 0.0) value.scale(1.0 / length);
    }
    public static void normalize(Point3f value) {
        float length = (float) Math.sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
        if (length > 0.0f) {
            value.x /= length;
            value.y /= length;
            value.z /= length;
        }
    }
    public static void normalizeSqr(Point3f value) {
        float lengthSquared = value.x * value.x + value.y * value.y + value.z * value.z;
        if (lengthSquared > 0.0f) {
            value.x /= lengthSquared;
            value.y /= lengthSquared;
            value.z /= lengthSquared;
        }
    }
    public static void normalizeSqr(Vector3d value) {
        double lengthSquared = value.x * value.x + value.y * value.y + value.z * value.z;
        if (lengthSquared > 0.0) {
            value.x /= lengthSquared;
            value.y /= lengthSquared;
            value.z /= lengthSquared;
        }
    }
    public static double length(Point3f value) {
        return Math.sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
    }
}
