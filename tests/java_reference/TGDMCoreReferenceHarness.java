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

/**
 * Selected MaCRUISE/VandyCRUISE TGDM calls on the synthetic base core maps.
 * This is a source-comparison harness, not a reimplementation of the JIST
 * application driver. The one-iteration wrapper values are recorded in the
 * fixture provenance.
 */
public final class TGDMCoreReferenceHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            throw new IllegalArgumentException("usage: TGDMCoreReferenceHarness INPUT_DIR OUTPUT_DIR");
        }
        File input = new File(args[0]);
        File output = new File(args[1]);
        if (!output.isDirectory() && !output.mkdirs()) {
            throw new IOException("cannot create output directory: " + output);
        }

        float[][][] membership = readNifti(new File(input, "initial_membership.nii.gz"));
        float[][][] wm = readNifti(new File(input, "white_matter.nii.gz"));
        float[][][] gm = readNifti(new File(input, "gray_matter.nii.gz"));
        float[][][] aceGm = readNifti(new File(input, "ace_gray_matter.nii.gz"));
        float[][][] skeleton = readNifti(new File(input, "ace_skeleton.nii.gz"));
        float[][][][] gvf = readNiftiVector(new File(input, "gvf.nii.gz"));
        int nx = membership.length, ny = membership[0].length, nz = membership[0][0].length;
        for (float[][][] values : new float[][][][] {wm, gm, aceGm, skeleton}) {
            if (values.length != nx || values[0].length != ny || values[0][0].length != nz) {
                throw new IllegalArgumentException("scalar input NIfTIs must have identical dimensions");
            }
        }
        if (gvf.length != nx || gvf[0].length != ny || gvf[0][0].length != nz) {
            throw new IllegalArgumentException("GVF dimensions must match scalar inputs");
        }

        ImageDataFloat initMembership = new ImageDataFloat(membership);
        float[][][] initial = initMembership.toArray3d();
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++)
                    initial[x][y][z] = 126.9f - initial[x][y][z];
        DistanceField distanceField = new DistanceField();
        ImageDataFloat signedSeed = distanceField.solve(new ImageDataFloat(initial), 8.0);
        MaGenericTGDM tgdm = new MaGenericTGDM(ConnectivityRule.CONNECT_18_6);
        ImageDataFloat smoothed = tgdm.solveSmoothSurface(signedSeed, 1.0f, 1);

        ImageDataFloat white = new ImageDataFloat(wm);
        ImageDataFloat gray = new ImageDataFloat(gm);
        ImageDataFloat aceGray = new ImageDataFloat(aceGm);
        ImageDataFloat thinSkeleton = new ImageDataFloat(skeleton);
        ImageDataFloat vector = new ImageDataFloat(gvf);

        ImageDataFloat firstInner = tgdm.solveInnerSurfaceFromOuter(
                smoothed, null, white, null, 0.2f, 0.0f, 1.0f, 126.9f, 0.6f, 1);

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
                        value = (wmP - 0.6f) / (1.0f - 0.6f);
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
                firstInner, null, new ImageDataFloat(sign), vector,
                0.15f, 1.0f, 1.5f, 1);
        ImageDataFloat centralShrink = tgdm.solveInnerSurfaceFromCentral(
                central, null, null, null, 0.15f, 0.0f, -1.0f,
                126.9f, 0.6f, 1);
        ImageDataFloat finalInner = tgdm.solveInnerSurfaceFromCentral(
                centralShrink, null, white, null, 0.2f, 0.0f, -1.0f,
                126.9f, 0.6f, 1);
        ImageDataFloat outer = tgdm.solveOuterSurface(
                central, null, gray, null, thinSkeleton,
                1.0f, 0.0f, 1.0f, 0.4f, 1);

        writeNifti(new File(output, "signed_seed.nii.gz"), signedSeed.toArray3d());
        writeNifti(new File(output, "smoothed_initial.nii.gz"), smoothed.toArray3d());
        writeNifti(new File(output, "initial_inner.nii.gz"), firstInner.toArray3d());
        writeNifti(new File(output, "tissue_sign.nii.gz"), sign);
        writeNifti(new File(output, "central.nii.gz"), central.toArray3d());
        writeNifti(new File(output, "central_shrink.nii.gz"), centralShrink.toArray3d());
        writeNifti(new File(output, "inner.nii.gz"), finalInner.toArray3d());
        writeNifti(new File(output, "outer.nii.gz"), outer.toArray3d());
        System.out.println("Wrote selected Java TGDM core reference outputs to " + output.getAbsolutePath());
    }

    static float[][][] readNifti(File file) throws IOException {
        try (InputStream input = new GZIPInputStream(new FileInputStream(file))) {
            byte[] headerBytes = input.readNBytes(352);
            if (headerBytes.length != 352) throw new IOException("short NIfTI header: " + file);
            ByteBuffer header = ByteBuffer.wrap(headerBytes).order(ByteOrder.LITTLE_ENDIAN);
            if (header.getInt(0) != 348 || header.getShort(40) != 3)
                throw new IOException("expected 3D NIfTI: " + file);
            int nx = header.getShort(42), ny = header.getShort(44), nz = header.getShort(46);
            int datatype = header.getShort(70);
            if (datatype != 2 && datatype != 8 && datatype != 16)
                throw new IOException("unsupported datatype: " + datatype);
            float[][][] values = new float[nx][ny][nz];
            byte[] voxel = new byte[datatype == 2 ? 1 : 4];
            for (int z = 0; z < nz; z++) for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
                int count = input.readNBytes(voxel, 0, voxel.length);
                if (count != voxel.length) throw new IOException("short NIfTI data: " + file);
                ByteBuffer decoded = ByteBuffer.wrap(voxel).order(ByteOrder.LITTLE_ENDIAN);
                values[x][y][z] = datatype == 2 ? (voxel[0] & 0xff)
                        : datatype == 8 ? decoded.getInt() : decoded.getFloat();
            }
            return values;
        }
    }

    /** NIfTI-1 vector intent with the 4th axis as three component planes. */
    static float[][][][] readNiftiVector(File file) throws IOException {
        try (InputStream input = new GZIPInputStream(new FileInputStream(file))) {
            byte[] headerBytes = input.readNBytes(352);
            if (headerBytes.length != 352) throw new IOException("short NIfTI header: " + file);
            ByteBuffer header = ByteBuffer.wrap(headerBytes).order(ByteOrder.LITTLE_ENDIAN);
            if (header.getInt(0) != 348 || header.getShort(40) != 4
                    || header.getShort(48) != 3 || header.getShort(70) != 16)
                throw new IOException("expected 4D float32 three-component NIfTI: " + file);
            int nx = header.getShort(42), ny = header.getShort(44), nz = header.getShort(46);
            float[][][][] values = new float[nx][ny][nz][3];
            byte[] voxel = new byte[4];
            for (int c = 0; c < 3; c++) for (int z = 0; z < nz; z++)
                for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
                    int count = input.readNBytes(voxel, 0, 4);
                    if (count != 4) throw new IOException("short vector NIfTI data: " + file);
                    values[x][y][z][c] = ByteBuffer.wrap(voxel).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                }
            return values;
        }
    }

    static void writeNifti(File file, float[][][] values) throws IOException {
        int nx = values.length, ny = values[0].length, nz = values[0][0].length;
        ByteBuffer header = ByteBuffer.allocate(352).order(ByteOrder.LITTLE_ENDIAN);
        header.putInt(0, 348); header.putShort(40, (short) 3);
        header.putShort(42, (short) nx); header.putShort(44, (short) ny); header.putShort(46, (short) nz);
        for (int index = 4; index < 8; index++) header.putShort(40 + 2 * index, (short) 1);
        header.putShort(70, (short) 16); header.putShort(72, (short) 32);
        header.putFloat(76, 1.0f); header.putFloat(80, 1.0f); header.putFloat(84, 1.0f); header.putFloat(88, 1.0f);
        header.putFloat(108, 352.0f); header.put(123, (byte) 2);
        header.putShort(252, (short) 1); header.putShort(254, (short) 1);
        header.putFloat(256, 1.0f); header.putFloat(280, 1.0f); header.putFloat(300, 1.0f); header.putFloat(320, 1.0f);
        header.put(344, (byte) 'n'); header.put(345, (byte) '+'); header.put(346, (byte) '1');
        try (GZIPOutputStream gzip = new GZIPOutputStream(new FileOutputStream(file))) {
            gzip.write(header.array());
            ByteBuffer voxel = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
            for (int z = 0; z < nz; z++) for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
                voxel.clear(); voxel.putFloat(values[x][y][z]); gzip.write(voxel.array());
            }
        }
    }
}
