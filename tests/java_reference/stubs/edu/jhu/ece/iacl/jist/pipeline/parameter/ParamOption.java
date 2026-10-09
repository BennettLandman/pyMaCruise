package edu.jhu.ece.iacl.jist.pipeline.parameter;
/** Headless option parameter stub. */
public class ParamOption {
    private int index;
    public ParamOption(String name, String[] options) {}
    public void setValue(int value) { index=value; }
    public int getIndex() { return index; }
    public void setMandatory(boolean value) {}
}
