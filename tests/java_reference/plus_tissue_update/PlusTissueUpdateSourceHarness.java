package edu.vanderbilt.masi.plugins.CRUISE.MaCRUISE;

import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

/** Runs the preserved selected CombineMAandTOADS source without the JIST UI. */
public final class PlusTissueUpdateSourceHarness {
    private static float[] readArray(DataInputStream in, int count) throws Exception {
        float[] values=new float[count]; for(int i=0;i<count;i++) values[i]=in.readFloat(); return values;
    }
    private static float[][][] toVolume(float[] flat,int nx,int ny,int nz) {
        float[][][] volume=new float[nx][ny][nz];
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++)
            volume[x][y][z]=flat[x+nx*y+nx*ny*z];
        return volume;
    }
    private static void put(CombineMAandTOADS algorithm,String name,Object value) throws Exception {
        Field field=CombineMAandTOADS.class.getDeclaredField(name); field.setAccessible(true); field.set(algorithm,value);
    }
    private static void write(Path path,ImageData image) throws Exception {
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(
                new GZIPOutputStream(new FileOutputStream(path.toFile()),1<<18)))) {
            out.writeInt(image.getRows());out.writeInt(image.getCols());out.writeInt(image.getSlices());
            for(int x=0;x<image.getRows();x++) for(int y=0;y<image.getCols();y++)
                for(int z=0;z<image.getSlices();z++) out.writeFloat(image.getFloat(x,y,z));
        }
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=2) throw new IllegalArgumentException("usage: INPUT_BIN OUTPUT_DIR");
        int nx,ny,nz; float[][][] gmMa,wmMa,brain,gmToads,wmToads,maskToads;
        try(DataInputStream in=new DataInputStream(new BufferedInputStream(new FileInputStream(args[0]),1<<20))) {
            nx=in.readInt();ny=in.readInt();nz=in.readInt();int count=nx*ny*nz;
            gmMa=toVolume(readArray(in,count),nx,ny,nz); wmMa=toVolume(readArray(in,count),nx,ny,nz);
            brain=toVolume(readArray(in,count),nx,ny,nz); gmToads=toVolume(readArray(in,count),nx,ny,nz);
            wmToads=toVolume(readArray(in,count),nx,ny,nz); maskToads=toVolume(readArray(in,count),nx,ny,nz);
            if(in.read()!=-1) throw new IllegalArgumentException("unexpected trailing bytes");
        }
        CombineMAandTOADS algorithm=new CombineMAandTOADS();
        put(algorithm,"GM_MA",volume(gmMa));put(algorithm,"WM_MA",volume(wmMa));put(algorithm,"BrainMask_MA",volume(brain));
        put(algorithm,"GM_toads",volume(gmToads));put(algorithm,"WM_toads",volume(wmToads));put(algorithm,"Mask_toads",volume(maskToads));
        ParamOption erosion=new ParamOption("Erode Method",new String[]{"3D","3Dmatlab"});erosion.setValue(0);put(algorithm,"ErodeMethod",erosion);
        ParamDouble confident=new ParamDouble("GMconfident level");confident.setValue(0.001);put(algorithm,"GM_confident_level",confident);
        ParamVolume gmOut=new ParamVolume("combined GM"),wmOut=new ParamVolume("combined WM"),maskOut=new ParamVolume("combined Mask");
        put(algorithm,"GM_combine",gmOut);put(algorithm,"WM_combine",wmOut);put(algorithm,"Mask_combine",maskOut);
        algorithm.execute(null);
        Path output=Path.of(args[1]);Files.createDirectories(output);
        write(output.resolve("java-gm.raw.gz"),gmOut.getImageData());
        write(output.resolve("java-wm.raw.gz"),wmOut.getImageData());
        write(output.resolve("java-mask.raw.gz"),maskOut.getImageData());
        System.out.println("shape="+nx+"x"+ny+"x"+nz+" source=CombineMAandTOADS.java erosion=3D GMconfident=0.001");
    }
    private static ParamVolume volume(float[][][] data) {
        ParamVolume parameter=new ParamVolume("source image");parameter.setValue(new ImageDataFloat(data));return parameter;
    }
}
