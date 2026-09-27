package edu.jhu.ece.iacl.jist.utility;

/** Logging is irrelevant to deterministic reference calculations. */
public final class JistLogger {
    public static final int INFO = 0;
    private JistLogger() {}
    public static void logOutput(int level, String message) {}
}
