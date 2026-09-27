package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Point3i;

/** Point type bridge required by the upstream connectivity helper classes. */
public final class GenericTGDM {
    private GenericTGDM() {}
    public static class NBPoint extends Point3i {}
}
