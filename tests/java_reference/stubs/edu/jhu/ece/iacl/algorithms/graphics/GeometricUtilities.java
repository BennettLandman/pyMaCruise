package edu.jhu.ece.iacl.algorithms.graphics;

import javax.vecmath.Vector3d;

/** Numerical helper subset used by the unmodified VandyCRUISE profiles. */
public final class GeometricUtilities {
    private GeometricUtilities() {}
    public static void normalize(Vector3d value) {
        double length = value.length();
        if (length > 0.0) value.scale(1.0 / length);
    }
}
