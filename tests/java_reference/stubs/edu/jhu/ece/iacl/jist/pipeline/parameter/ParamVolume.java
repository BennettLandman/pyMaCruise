package edu.jhu.ece.iacl.jist.pipeline.parameter;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
/** Headless image parameter stub. */
public class ParamVolume {
    private ImageData value;
    public ParamVolume(String name) {}
    public ImageData getImageData() { return value; }
    public void setValue(ImageData image) { value=image; }
}
