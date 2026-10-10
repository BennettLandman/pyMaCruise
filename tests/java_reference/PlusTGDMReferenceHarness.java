import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.vanderbilt.masi.plugins.CRUISE.tgdm.MaGenericTGDM;

/**
 * Selected MaGenericTGDM plus-wrapper call sequence on held-fixed Python
 * inputs. This compares the nested surface core only; it does not reimplement
 * the MaCRUISEplus workflow or its upstream TOADS/MaSeg stages.
 */
public final class PlusTGDMReferenceHarness {
    private static float[][][] readVolume(DataInputStream in, int nx, int ny, int nz)
            throws Exception {
        float[][][] a = new float[nx][ny][nz];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
            a[x][y][z] = in.readFloat();
        return a;
    }

    private static void writeVolume(Path path, float[][][] a) throws Exception {
        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(
                new GZIPOutputStream(new FileOutputStream(path.toFile()), 1 << 18)))) {
            for (int x=0;x<a.length;x++) for (int y=0;y<a[0].length;y++)
                for (int z=0;z<a[0][0].length;z++) out.writeFloat(a[x][y][z]);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2 || args.length > 3) throw new IllegalArgumentException("usage: INPUT_BIN OUTPUT_DIR [precomputed_initial_level_set]");
        boolean precomputedInitial = args.length == 3 && args[2].equals("precomputed_initial_level_set");
        if (args.length == 3 && !precomputedInitial) throw new IllegalArgumentException("unknown input mode: " + args[2]);
        Path output = Path.of(args[1]); Files.createDirectories(output);
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(
                new FileInputStream(args[0]), 1 << 20))) {
            int nx=in.readInt(), ny=in.readInt(), nz=in.readInt();
            float[][][] membership=readVolume(in,nx,ny,nz);
            float[][][] wm=readVolume(in,nx,ny,nz);
            float[][][] gm=readVolume(in,nx,ny,nz);
            float[][][] aceGm=readVolume(in,nx,ny,nz);
            float[][][] skeleton=readVolume(in,nx,ny,nz);
            float[][][][] gvf=new float[nx][ny][nz][3];
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                for (int c=0;c<3;c++) gvf[x][y][z][c]=in.readFloat();
            if (in.read() != -1) throw new IllegalArgumentException("unexpected trailing bytes");

            ImageDataFloat initialLevel;
            if (precomputedInitial) {
                initialLevel = new ImageDataFloat(membership);
            } else {
                ImageDataFloat initialMembership=new ImageDataFloat(membership);
                float[][][] initial=initialMembership.toArray3d();
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                    initial[x][y][z]=126.9f-initial[x][y][z];
                initialLevel=new DistanceField().solve(new ImageDataFloat(initial),8.0);
            }
            MaGenericTGDM tgdm=new MaGenericTGDM(ConnectivityRule.CONNECT_18_6);
            ImageDataFloat smoothed=tgdm.solveSmoothSurface(initialLevel,1.0f,6);
            ImageDataFloat white=new ImageDataFloat(wm), gray=new ImageDataFloat(gm);
            ImageDataFloat ace=new ImageDataFloat(aceGm), thin=new ImageDataFloat(skeleton);
            ImageDataFloat vector=new ImageDataFloat(gvf);

            ImageDataFloat firstInner=tgdm.solveInnerSurfaceFromOuter(
                    smoothed,null,white,null,0.2f,0.0f,1.0f,126.9f,0.5f,4);
            float[][][] sign=new float[nx][ny][nz];
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                float wmP=wm[x][y][z]/255.0f, gmP=aceGm[x][y][z]/255.0f;
                float csfP=1.0f-wmP-gmP, value;
                if (wmP>csfP) { value=(wmP-0.5f)/(1.0f-0.5f); if(value<0) value=0; }
                else { value=(1.0f-csfP-0.4f)/0.4f; if(value>0) value=0; }
                int q=(int)(value*127+127); sign[x][y][z]=q<254?q:254;
            }
            ImageDataFloat central=tgdm.solveCentralSurface(
                    firstInner,null,new ImageDataFloat(sign),vector,0.15f,1.0f,1.5f,7);
            ImageDataFloat centralShrink=tgdm.solveInnerSurfaceFromCentral(
                    central,null,null,null,0.15f,0.0f,-1.0f,126.9f,0.5f,1);
            ImageDataFloat inner=tgdm.solveInnerSurfaceFromCentral(
                    centralShrink,null,white,null,0.2f,0.0f,-1.0f,126.9f,0.5f,4);
            ImageDataFloat outer=tgdm.solveOuterSurface(
                    central,null,gray,null,thin,1.0f,0.0f,1.0f,0.4f,7);
            writeVolume(output.resolve("initial_level_set.raw.gz"),initialLevel.toArray3d());
            writeVolume(output.resolve("initial_inner.raw.gz"),firstInner.toArray3d());
            writeVolume(output.resolve("central.raw.gz"),central.toArray3d());
            writeVolume(output.resolve("inner.raw.gz"),inner.toArray3d());
            writeVolume(output.resolve("outer.raw.gz"),outer.toArray3d());
        }
    }
}
