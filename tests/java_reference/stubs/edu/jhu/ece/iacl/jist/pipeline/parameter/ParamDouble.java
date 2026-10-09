package edu.jhu.ece.iacl.jist.pipeline.parameter;
/** Headless numeric parameter stub. */
public class ParamDouble {
    private double value;
    public ParamDouble(String name) {}
    public void setValue(double value) { this.value=value; }
    public float getFloat() { return (float)value; }
    public void setMandatory(boolean mandatory) {}
}
