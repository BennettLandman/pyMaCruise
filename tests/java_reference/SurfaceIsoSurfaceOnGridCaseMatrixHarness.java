import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

import java.io.FileWriter;

/**
 * Development-only case oracle for the selected JIST isosurface implementation.
 *
 * <p>Runs all 256 binary corner configurations through the supplied
 * {@code IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)} body.
 * The source is authored by Blake Lucas; see the source inventory and
 * attribution notes in {@code license_info.md}. This harness and its JIST
 * stubs are never imported by the Python runtime.</p>
 */
public final class SurfaceIsoSurfaceOnGridCaseMatrixHarness {
    private static void writeArray(FileWriter writer, float[][] values) throws Exception {
        writer.write("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) writer.write(",");
            writer.write("[" + values[i][0] + "," + values[i][1] + "," + values[i][2] + "]");
        }
        writer.write("]");
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 1 || args.length > 2) throw new IllegalArgumentException(
            "usage: SurfaceIsoSurfaceOnGridCaseMatrixHarness output.json [variants]");
        int variantCount = args.length == 2 ? Integer.parseInt(args[1]) : 1;
        if (variantCount < 1 || variantCount > 16)
            throw new IllegalArgumentException("variants must be between 1 and 16");
        try (FileWriter writer = new FileWriter(args[0])) {
            writer.write("{\"source_method\":\"solveOriginal(ImageDataFloat,18,0,false)\",\"patterns\":[");
            for (int pattern = 0; pattern < 256; pattern++) {
                if (pattern > 0) writer.write(",");
                writer.write("{\"pattern\":" + pattern + ",\"variants\":[");
                for (int variant = 0; variant < variantCount; variant++) {
                    if (variant > 0) writer.write(",");
                    float[][][] values = new float[2][2][2];
                    float[] cornerValues = new float[8];
                    for (int z = 0; z < 2; z++) {
                        for (int y = 0; y < 2; y++) {
                            for (int x = 0; x < 2; x++) {
                                int bit = (z << 2) | (y << 1) | x;
                                float magnitude = 0.31f + (variant + 1) * 0.47f
                                    + (bit + 1) * (0.19f + variant * 0.07f)
                                    + (bit + 1) * (bit + 1) * (0.011f + variant * 0.003f);
                                values[x][y][z] = ((pattern >> bit) & 1) == 0
                                    ? -magnitude : magnitude;
                                cornerValues[bit] = values[x][y][z];
                            }
                        }
                    }
                    EmbeddedSurface surface = new IsoSurfaceOnGrid().solveOriginal(
                        new ImageDataFloat(values), 18, 0.0f, false);
                    writer.write("{\"corner_values\":[");
                    for (int i = 0; i < cornerValues.length; i++) {
                        if (i > 0) writer.write(",");
                        writer.write(Float.toString(cornerValues[i]));
                    }
                    writer.write("],\"vertices\":");
                    if (surface == null) {
                        writer.write("[],\"faces\":[]}");
                        continue;
                    }
                    float[][] vertices = new float[surface.points.length][3];
                    for (int i = 0; i < surface.points.length; i++) {
                        vertices[i][0] = surface.points[i].x;
                        vertices[i][1] = surface.points[i].y;
                        vertices[i][2] = surface.points[i].z;
                    }
                    writeArray(writer, vertices);
                    writer.write(",\"faces\":[");
                    for (int i = 0; i < surface.indices.length; i += 3) {
                        if (i > 0) writer.write(",");
                        writer.write("[" + surface.indices[i] + "," + surface.indices[i + 1]
                            + "," + surface.indices[i + 2] + "]");
                    }
                    writer.write("]}");
                }
                writer.write("]}");
            }
            writer.write("]}\n");
        }
    }
}
