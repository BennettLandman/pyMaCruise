package gov.nih.mipav.model.algorithms;

import java.util.BitSet;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.view.ViewJProgressBar;

/** Minimal headless API shim used only to compile the selected reference algorithm. */
public class AlgorithmBase {
    protected ModelImage srcImage;
    protected BitSet mask = new BitSet();
    protected boolean threadStopped;
    private boolean completed;

    public AlgorithmBase(Object parent, ModelImage source) { srcImage = source; }
    public void finalize() { }
    protected void setStartTime() { }
    protected void computeElapsedTime() { }
    protected void fireProgressStateChanged(String message) { }
    protected void fireProgressStateChanged(String a, String b) { }
    protected void displayError(String message) { }
    protected ViewJProgressBar getProgressChangeListener() { return null; }
    protected void errorCleanUp(String message, boolean flag) { }
    protected void setCompleted(boolean value) { completed = value; }
    public boolean isCompleted() { return completed; }
}
