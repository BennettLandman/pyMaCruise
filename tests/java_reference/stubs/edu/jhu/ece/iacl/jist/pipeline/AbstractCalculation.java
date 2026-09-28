package edu.jhu.ece.iacl.jist.pipeline;

/** Minimal test-only shell for the numerical GVF source dependencies. */
public class AbstractCalculation {
    public AbstractCalculation() {}
    public AbstractCalculation(AbstractCalculation parent) {}
    public void setLabel(String label) {}
    public void setTotalUnits(int units) {}
    public void setTotalUnits(long units) {}
    public void incrementCompletedUnits() {}
    public void decrementCompletedUnits() {}
    public void markCompleted() {}
}
