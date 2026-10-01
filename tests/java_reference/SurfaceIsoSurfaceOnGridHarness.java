import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** Run the selected IsoSurfaceOnGrid Java body on a little-endian float volume. */
public final class SurfaceIsoSurfaceOnGridHarness {
    private static String sha256(byte[] bytes) throws Exception {
        byte[] digest = MessageDigest.getInstance("SHA-256").digest(bytes);
        StringBuilder out = new StringBuilder();
        for (byte value : digest) out.append(String.format("%02x", value & 0xff));
        return out.toString();
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 5) throw new IllegalArgumentException(
            "usage: SurfaceIsoSurfaceOnGridHarness raw-f32 nx ny nz output.json");
        Path input = Path.of(args[0]);
        int nx = Integer.parseInt(args[1]);
        int ny = Integer.parseInt(args[2]);
        int nz = Integer.parseInt(args[3]);
        Path output = Path.of(args[4]);
        byte[] bytes = Files.readAllBytes(input);
        if (bytes.length != nx * ny * nz * 4) throw new IllegalArgumentException("raw input byte count mismatch");
        ByteBuffer values = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[][][] volume = new float[nx][ny][nz];
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++) volume[x][y][z] = values.getFloat();

        ImageDataFloat image = new ImageDataFloat(volume);
        EmbeddedSurface surface = new IsoSurfaceOnGrid().solveOriginal(image, 18, 0.0f, false);
        if (surface == null) throw new IllegalStateException("selected extractor returned no surface");
        try (FileWriter writer = new FileWriter(output.toFile())) {
            writer.write("{\"input_sha256\":\"" + sha256(bytes) + "\",\"dimensions\":[" + nx + "," + ny + "," + nz + "],\"connectivity\":18,\"iso_level\":0,\"inclusive\":false,\"vertices\":[");
            for (int i = 0; i < surface.points.length; i++) {
                if (i > 0) writer.write(",");
                writer.write("[" + Float.toString(surface.points[i].x) + "," + Float.toString(surface.points[i].y) + "," + Float.toString(surface.points[i].z) + "]");
            }
            writer.write("],\"faces\":[");
            for (int i = 0; i < surface.indices.length; i += 3) {
                if (i > 0) writer.write(",");
                writer.write("[" + surface.indices[i] + "," + surface.indices[i + 1] + "," + surface.indices[i + 2] + "]");
            }
            writer.write("]}\n");
        }
    }
}
