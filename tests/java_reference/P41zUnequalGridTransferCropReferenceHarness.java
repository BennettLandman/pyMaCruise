import edu.jhmi.rad.medic.utilities.CropParameters;
import edu.jhmi.rad.medic.utilities.CubicVolumeCropper;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.ViewUserInterface;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;

/** Replay CropBasedOnGM on a full GM reference and smaller transferred targets. */
public final class P41zUnequalGridTransferCropReferenceHarness {
    private static ImageData read(Path path, int x, int y, int z, boolean integer) throws Exception {
        ImageData image = integer ? new ImageDataInt(x, y, z) : new ImageDataFloat(x, y, z);
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(Files.newInputStream(path)))) {
            for (int i=0;i<x;i++) for (int j=0;j<y;j++) for (int k=0;k<z;k++) {
                if (integer) image.set(i,j,k,in.readInt()); else image.set(i,j,k,in.readFloat());
            }
            if (in.read()!=-1) throw new IllegalArgumentException("input contains trailing data: "+path);
        }
        return image;
    }
    private static void write(Path path, ImageData image) throws Exception {
        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(path)))) {
            for (int i=0;i<image.getRows();i++) for (int j=0;j<image.getCols();j++) for (int k=0;k<image.getSlices();k++) out.writeFloat(image.getFloat(i,j,k));
        }
    }
    public static void main(String[] args) throws Exception {
        if (args.length != 9) throw new IllegalArgumentException("usage: <input-dir> <output-dir> <gm-x> <gm-y> <gm-z> <target-x> <target-y> <target-z> <edge>");
        Path input=Path.of(args[0]), output=Path.of(args[1]); Files.createDirectories(output);
        int gx=Integer.parseInt(args[2]),gy=Integer.parseInt(args[3]),gz=Integer.parseInt(args[4]);
        int tx=Integer.parseInt(args[5]),ty=Integer.parseInt(args[6]),tz=Integer.parseInt(args[7]),edge=Integer.parseInt(args[8]);
        ViewUserInterface.create(); Path library=Path.of(System.getProperty("user.home"),"jist-library"); Files.createDirectories(library); Preferences.setProperty("JIST_Library_Directory",library.toString());
        CubicVolumeCropper cropper=new CubicVolumeCropper();
        ImageData gm=read(input.resolve("atlas_gm.i32be"),gx,gy,gz,true); ImageData croppedGm=null; CropParameters parameters;
        try { croppedGm=cropper.crop(gm,0.0,edge); parameters=cropper.getLastCropParams(); write(output.resolve("atlas_gm.f32be"),croppedGm); }
        finally { if(croppedGm!=null)croppedGm.dispose(); gm.dispose(); }
        String[] names={"csf","gm","wm_fill","wm_mask"};
        for(String name:names){
            ImageData target=read(input.resolve(name+".f32be"),tx,ty,tz,false); ImageData cropped=null;
            try { cropped=cropper.crop(target,parameters); write(output.resolve(name+".f32be"),cropped); }
            finally { if(cropped!=null)cropped.dispose(); target.dispose(); }
        }
        Files.writeString(output.resolve("bounds.txt"), parameters.xmin+","+parameters.xmax+";"+parameters.ymin+","+parameters.ymax+";"+parameters.zmin+","+parameters.zmax+";"+parameters.getCroppedRows()+","+parameters.getCroppedCols()+","+parameters.getCroppedSlices()+"\n");
    }
}
