package gov.nih.mipav.view;

import java.util.StringTokenizer;

/** Numeric token parsing used by the headless registration reference harness. */
public final class MipavUtil {
    private MipavUtil() {}
    public static int getInt(StringTokenizer value) { return Integer.parseInt(value.nextToken().trim()); }
    public static float getFloat(StringTokenizer value) { return Float.parseFloat(value.nextToken().trim()); }
    public static boolean getBoolean(StringTokenizer value) { return Boolean.parseBoolean(value.nextToken().trim()); }
}
