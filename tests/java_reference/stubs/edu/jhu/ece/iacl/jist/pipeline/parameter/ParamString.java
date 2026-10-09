package edu.jhu.ece.iacl.jist.pipeline.parameter;
/** Headless string parameter stub required by selected JIST source. */
public class ParamString {
    private String value;
    public ParamString() {}
    public ParamString(String name, String initialValue) { value=initialValue; }
    public String getValue() { return value; }
    public void setValue(String value) { this.value=value; }
}
