package gov.nih.mipav.view;

/** Headless progress-listener shim; reference harness runs with a null listener. */
public class ViewJProgressBar {
    public void updateValue(int value, boolean separateThread) { }
    public void updateValue(int value) { }
}
