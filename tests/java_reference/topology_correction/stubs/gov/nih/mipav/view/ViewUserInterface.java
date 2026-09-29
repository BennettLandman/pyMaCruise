package gov.nih.mipav.view;

/** Headless UI singleton shim. */
public class ViewUserInterface {
    private static final ViewUserInterface INSTANCE = new ViewUserInterface();
    public static ViewUserInterface getReference() { return INSTANCE; }
}
