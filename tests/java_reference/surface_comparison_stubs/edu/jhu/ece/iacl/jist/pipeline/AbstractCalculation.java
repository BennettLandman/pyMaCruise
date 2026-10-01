package edu.jhu.ece.iacl.jist.pipeline;

/** Numerical harness surface for the progress-only superclass API. */
public abstract class AbstractCalculation {
    protected void setLabel(String value) {}
    protected void setTotalUnits(int value) {}
    protected void incrementCompletedUnits() {}
    protected void markCompleted() {}
}
