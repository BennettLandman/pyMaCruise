package edu.jhu.ece.iacl.algorithms.volume;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.vanderbilt.masi.algorithms.CRUISE.ace.DistanceFieldyh;

/**
 * Validation-only bridge for the omitted DistanceField class.
 * DistanceFieldyh is supplied with the source tree and documents that it is
 * identical to DistanceField; no production Python dependency uses this.
 */
public final class DistanceField extends AbstractCalculation {
    public DistanceField() {
        super();
    }

    public DistanceField(AbstractCalculation parent) {
        super(parent);
    }

    public ImageDataFloat solve(ImageDataFloat volume, double maxDistance) {
        return DistanceFieldyh.doSolve(volume, maxDistance);
    }
}
