package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

/** Compile-time bridge for FastMarchingGradient's unused list overload. */
public final class TrackingTGDM {
    private TrackingTGDM() {}
    public static class EmbeddedNBPoint extends GenericTGDM.NBPoint {
        public Vector3d vec;
    }
}
