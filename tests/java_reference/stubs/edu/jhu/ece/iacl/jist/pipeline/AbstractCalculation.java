package edu.jhu.ece.iacl.jist.pipeline;

/** Minimal non-GUI runner API needed to compile the original TGDM source. */
public class AbstractCalculation {
    public AbstractCalculation() {}
    public AbstractCalculation(AbstractCalculation parent) {}
    public void setLabel(String value) {}
    public void setTotalUnits(int value) {}
    public void incrementCompletedUnits() {}
    public void decrementCompletedUnits() {}
    public void markCompleted() {}
}
