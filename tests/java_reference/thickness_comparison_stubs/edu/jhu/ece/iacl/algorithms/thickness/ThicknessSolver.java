package edu.jhu.ece.iacl.algorithms.thickness;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;

/** Minimal superclass state/API shell; selected subclass methods compile unchanged. */
public abstract class ThicknessSolver extends AbstractCalculation {
    protected double isoVal = 0;
    protected double lambda = 0.25;
    protected double lagrangeDeltaT = 0.1;
    protected boolean ADJUST_BOUNDARY = false;
    protected int rows, cols, slices;
    protected EmbeddedSurface innerMesh, centralMesh, outerMesh;
    protected ImageData innerVol, outerVol;

    public EmbeddedSurface getInnerSurface() { return innerMesh; }
    public EmbeddedSurface getCentralSurface() { return centralMesh; }
    public EmbeddedSurface getOuterSurface() { return outerMesh; }
}
