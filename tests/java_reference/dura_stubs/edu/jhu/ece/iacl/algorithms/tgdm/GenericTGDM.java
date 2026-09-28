package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

/**
 * State-only superclass seam for the one-point DynamicPressureTGDM reference.
 * The selected GenericTGDM evolution loop is deliberately not represented here.
 */
public class GenericTGDM extends AbstractCalculation {
    public static class NBPoint extends Point3i {
        public double w;
        public int index;
    }

    protected int rows, cols, slices;
    protected float[][][] innerBound, outerBound, regionForces, skeletalBound, implicitLevelSet;
    protected float[][][][] advectionForces;
    protected Vector3d curvatureGrad = new Vector3d();
    protected Vector3d regionGrad = new Vector3d();
    protected Vector3d advectionGrad = new Vector3d();
    protected Vector3d advectionForce = new Vector3d();
    protected DistanceField df;
    protected int connectivityRule;
    protected double timeStep;

    public GenericTGDM(int rule) { this.connectivityRule = rule; }
    public GenericTGDM(AbstractCalculation parent, int rule) {
        super(parent);
        this.connectivityRule = rule;
    }
    protected void execute(ProfileTGDM profile) {}
}
