import java.io.File;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;

import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.vanderbilt.masi.plugins.CRUISE.tgdm.MaGenericTGDM;

/** Run the supplied MaCRUISE/VandyCRUISE solver on the deterministic phantom. */
public final class TGDMReferenceHarness {
    private static final int SIZE = 25;
    private static final float CENTER = 12.0f;

    public static void main(String[] args) throws Exception {
        if (args.length < 1 || args.length > 3) {
            throw new IllegalArgumentException("usage: TGDMReferenceHarness OUTPUT_DIR [INPUT_DIR [WM_ISO]]");
        }
        float wmIso = args.length > 2 ? Float.parseFloat(args[2]) : 0.5f;
        File output = new File(args[0]);
        if (!output.isDirectory() && !output.mkdirs()) {
            throw new IOException("cannot create output directory: " + output);
        }

        float[][][] membership;
        float[][][] wm;
        float[][][] gm;
        float[][][] aceGm;
        float[][][] skeleton;
        float[][][] brainBound = null;
        int nx = SIZE, ny = SIZE, nz = SIZE;
        if (args.length >= 2) {
            File input = new File(args[1]);
            membership = readNifti(new File(input, "initial_membership.nii.gz"));
            wm = readNifti(new File(input, "white_matter.nii.gz"));
            gm = readNifti(new File(input, "gray_matter.nii.gz"));
            aceGm = readNifti(new File(input, "ace_gray_matter.nii.gz"));
            skeleton = readNifti(new File(input, "ace_skeleton.nii.gz"));
            File boundFile = new File(input, "tgdm_brain_bound.nii.gz");
            if (boundFile.isFile()) brainBound = readNifti(boundFile);
            nx = membership.length;
            ny = membership[0].length;
            nz = membership[0][0].length;
            for (float[][][] values : new float[][][][] {wm, gm, aceGm, skeleton}) {
                if (values.length != nx || values[0].length != ny || values[0][0].length != nz) {
                    throw new IllegalArgumentException("input NIfTIs must have identical dimensions");
                }
            }
        } else {
            membership = new float[SIZE][SIZE][SIZE];
            wm = new float[SIZE][SIZE][SIZE];
            gm = new float[SIZE][SIZE][SIZE];
            aceGm = new float[SIZE][SIZE][SIZE];
            skeleton = new float[SIZE][SIZE][SIZE];
        }
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    if (args.length >= 2) continue;
                    double dx = x - CENTER;
                    double dy = y - CENTER;
                    double dz = z - CENTER;
                    double radius = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    membership[x][y][z] = radius <= 8 ? 255 : 0;
                    wm[x][y][z] = radius <= 5 ? 255 : 0;
                    gm[x][y][z] = radius > 5 && radius <= 8 ? 220 : 0;
                    aceGm[x][y][z] = gm[x][y][z];
                    skeleton[x][y][z] = Math.abs(radius - 6.5) < 0.5 && z == 12 ? 254 : 0;
                }
            }
        }

        ImageDataFloat initMembership = new ImageDataFloat(membership);
        float[][][] initial = initMembership.toArray3d();
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    initial[x][y][z] = 126.9f - initial[x][y][z];
                }
            }
        }
        DistanceField distanceField = new DistanceField();
        ImageDataFloat signedSeed = distanceField.solve(new ImageDataFloat(initial), 8.0);
        MaGenericTGDM tgdm = new MaGenericTGDM(ConnectivityRule.CONNECT_18_6);
        ImageDataFloat smoothed = tgdm.solveSmoothSurface(signedSeed, 1.0f, 6);

        writeNifti(new File(output, "signed_seed.nii.gz"), signedSeed.toArray3d());
        writeNifti(new File(output, "smoothed_initial.nii.gz"), smoothed.toArray3d());

        ImageDataFloat white = new ImageDataFloat(wm);
        ImageDataFloat gray = new ImageDataFloat(gm);
        ImageDataFloat aceGray = new ImageDataFloat(aceGm);
        ImageDataFloat thinSkeleton = new ImageDataFloat(skeleton);
        ImageDataFloat gvf = new ImageDataFloat(new float[nx][ny][nz][3]);
        ImageDataFloat firstInner = tgdm.solveInnerSurfaceFromOuter(
                smoothed, brainBound == null ? null : new ImageDataFloat(brainBound),
                white, null, 0.2f, 0.0f, 1.0f, 126.9f, wmIso, 4);

        float[][][] sign = new float[nx][ny][nz];
        float[][][] wmValues = white.toArray3d();
        float[][][] aceValues = aceGray.toArray3d();
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    float wmP = wmValues[x][y][z] / 255.0f;
                    float gmP = aceValues[x][y][z] / 255.0f;
                    float csfP = 1.0f - wmP - gmP;
                    float value;
                    if (wmP > csfP) {
                        value = (wmP - 0.5f) / (1.0f - 0.5f);
                        if (value < 0) value = 0;
                    } else {
                        value = (1.0f - csfP - 0.4f) / 0.4f;
                        if (value > 0) value = 0;
                    }
                    int quantized = (int) (value * 127 + 127);
                    sign[x][y][z] = quantized < 254 ? quantized : 254;
                }
            }
        }

        ImageDataFloat central = tgdm.solveCentralSurface(
                firstInner, brainBound == null ? null : new ImageDataFloat(brainBound),
                new ImageDataFloat(sign), gvf,
                0.15f, 1.0f, 1.5f, 7);
        ImageDataFloat centralShrink = tgdm.solveInnerSurfaceFromCentral(
                central, null, null, null, 0.15f, 0.0f, -1.0f,
                126.9f, wmIso, 1);
        ImageDataFloat finalInner = tgdm.solveInnerSurfaceFromCentral(
                centralShrink, null, white, null, 0.2f, 0.0f, -1.0f,
                126.9f, wmIso, 4);
        ImageDataFloat outer = tgdm.solveOuterSurface(
                central, brainBound == null ? null : new ImageDataFloat(brainBound),
                gray, null, thinSkeleton,
                1.0f, 0.0f, 1.0f, 0.4f, 7);

        writeNifti(new File(output, "initial_inner.nii.gz"), firstInner.toArray3d());
        writeNifti(new File(output, "tissue_sign.nii.gz"), sign);
        writeNifti(new File(output, "central.nii.gz"), central.toArray3d());
        writeNifti(new File(output, "central_shrink.nii.gz"), centralShrink.toArray3d());
        writeNifti(new File(output, "inner.nii.gz"), finalInner.toArray3d());
        writeNifti(new File(output, "outer.nii.gz"), outer.toArray3d());
        System.out.println("Wrote TGDM Java reference outputs to " + output.getAbsolutePath());
    }

    static void writeNifti(File file, float[][][] values) throws IOException {
        int nx = values.length, ny = values[0].length, nz = values[0][0].length;
        ByteBuffer header = ByteBuffer.allocate(352).order(ByteOrder.LITTLE_ENDIAN);
        header.putInt(0, 348);
        header.putShort(40, (short) 3);
        header.putShort(42, (short) nx);
        header.putShort(44, (short) ny);
        header.putShort(46, (short) nz);
        for (int index = 4; index < 8; index++) header.putShort(40 + 2 * index, (short) 1);
        header.putShort(70, (short) 16); // NIFTI_TYPE_FLOAT32
        header.putShort(72, (short) 32);
        header.putFloat(76, 1.0f);
        header.putFloat(80, 1.0f);
        header.putFloat(84, 1.0f);
        header.putFloat(88, 1.0f);
        header.putFloat(108, 352.0f);
        header.put(123, (byte) 2); // millimetres
        header.putShort(252, (short) 1);
        header.putShort(254, (short) 1);
        header.putFloat(256, 1.0f);
        header.putFloat(280, 1.0f);
        header.putFloat(300, 1.0f);
        header.putFloat(320, 1.0f);
        header.put(344, (byte) 'n');
        header.put(345, (byte) '+');
        header.put(346, (byte) '1');

        try (GZIPOutputStream gzip = new GZIPOutputStream(new FileOutputStream(file))) {
            gzip.write(header.array());
            byte[] floatBytes = new byte[4];
            for (int z = 0; z < nz; z++) {
                for (int y = 0; y < ny; y++) {
                    for (int x = 0; x < nx; x++) {
                        int bits = Float.floatToIntBits(values[x][y][z]);
                        floatBytes[0] = (byte) (bits & 0xff);
                        floatBytes[1] = (byte) ((bits >>> 8) & 0xff);
                        floatBytes[2] = (byte) ((bits >>> 16) & 0xff);
                        floatBytes[3] = (byte) ((bits >>> 24) & 0xff);
                        gzip.write(floatBytes);
                    }
                }
            }
        }
    }

    /** Read little-endian 3D uint8/float32 RAS NIfTIs used by synthetic fixtures. */
    static float[][][] readNifti(File file) throws IOException {
        try (InputStream input = new GZIPInputStream(new FileInputStream(file))) {
            byte[] headerBytes = input.readNBytes(352);
            if (headerBytes.length != 352) throw new IOException("short NIfTI header: " + file);
            ByteBuffer header = ByteBuffer.wrap(headerBytes).order(ByteOrder.LITTLE_ENDIAN);
            if (header.getInt(0) != 348 || header.getShort(40) != 3) {
                throw new IOException("expected little-endian 3D NIfTI: " + file);
            }
            int nx = header.getShort(42), ny = header.getShort(44), nz = header.getShort(46);
            int datatype = header.getShort(70);
            if (datatype != 2 && datatype != 16) throw new IOException("unsupported NIfTI datatype " + datatype);
            float[][][] values = new float[nx][ny][nz];
            byte[] voxel = new byte[datatype == 2 ? 1 : 4];
            for (int z = 0; z < nz; z++) {
                for (int y = 0; y < ny; y++) {
                    for (int x = 0; x < nx; x++) {
                        int count = input.readNBytes(voxel, 0, voxel.length);
                        if (count != voxel.length) throw new IOException("short NIfTI data: " + file);
                        values[x][y][z] = datatype == 2
                            ? (voxel[0] & 0xff)
                            : ByteBuffer.wrap(voxel).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                    }
                }
            }
            return values;
        }
    }
}
