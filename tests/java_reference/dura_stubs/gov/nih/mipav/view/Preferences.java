package gov.nih.mipav.view;

import java.io.File;

/** Headless preference shim for loading the original JIST image algorithms. */
public final class Preferences {
    public static final String PREF_MAX_HEAP_SIZE = "PREF_MAX_HEAP_SIZE";
    private static final File LIBRARY = new File("build/java-reference/jist-empty-lib");

    public static String getProperty(String key) {
        if ("JIST_Library_Directory".equals(key)) {
            LIBRARY.mkdirs();
            return LIBRARY.getAbsolutePath();
        }
        return null;
    }

    public static void setProperty(String key, String value) { }
}
