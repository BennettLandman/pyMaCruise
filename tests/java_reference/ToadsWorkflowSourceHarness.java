import edu.jhmi.rad.medic.algorithms.AlgorithmMultipleToads;
import edu.jhmi.rad.medic.methods.MultipleToadDeformableAtlas;
import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.ModelStorageBase;
import java.io.*;
import java.nio.*;
import java.nio.file.*;
import java.security.*;
import java.util.*;

/**
 * Runs the selected V3 TOADS wrapper settings against the supplied source
 * atlas and a held-fixed float32 T1 volume. This test-only adapter deliberately
 * exposes only the terminal classification and its source-order hash.
 */
public final class ToadsWorkflowSourceHarness {
    private static String sha256(byte[] values) throws Exception {
        byte[] digest=MessageDigest.getInstance("SHA-256").digest(values);
        StringBuilder text=new StringBuilder();
        for(byte value:digest)text.append(String.format("%02x",value&255));
        return text.toString();
    }

    public static void main(String[] args) throws Exception {
        if(args.length!=10)throw new IllegalArgumentException(
                "input.raw atlas.txt nx ny nz sx sy sz maxIterations output.raw");
        int nx=Integer.parseInt(args[2]),ny=Integer.parseInt(args[3]),nz=Integer.parseInt(args[4]);
        int maxIterations=Integer.parseInt(args[8]);
        byte[] raw=Files.readAllBytes(Paths.get(args[0]));
        if(raw.length!=4*nx*ny*nz)throw new IllegalArgumentException("wrong raw T1 length");
        ByteBuffer buffer=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        float[] values=new float[nx*ny*nz];
        for(int i=0;i<values.length;i++)values[i]=buffer.getFloat();
        ModelImage image=new ModelImage(ModelStorageBase.FLOAT,new int[]{nx,ny,nz},"synthetic T1");
        image.importData(0,values,true);
        FileInfoBase info=image.getFileInfo(0);
        info.setResolutions(new float[]{Float.parseFloat(args[5]),Float.parseFloat(args[6]),Float.parseFloat(args[7])});
        info.setUnitsOfMeasure(new int[]{FileInfoBase.MILLIMETERS,FileInfoBase.MILLIMETERS,FileInfoBase.MILLIMETERS});
        info.setAxisOrientation(new int[]{FileInfoBase.ORI_L2R_TYPE,FileInfoBase.ORI_P2A_TYPE,FileInfoBase.ORI_I2S_TYPE});
        info.setImageOrientation(FileInfoBase.AXIAL);
        image.setImageOrientation(FileInfoBase.AXIAL);

        MultipleToadDeformableAtlas atlas=new MultipleToadDeformableAtlas(args[1]);
        AlgorithmMultipleToads algorithm=new AlgorithmMultipleToads(
                new ModelImage[]{image},1,new String[]{"T1_MPRAGE"},args[1],atlas,"hard_segmentation",
                0.2f,maxIterations,1.0e-4f,3.0f,0.5f,2.0f,0.2f,"rigid",
                4,50,5,0,true,"Chebyshev",3,30.0f,"18/6");
        algorithm.runAlgorithm();
        ModelImage[] results=algorithm.getResultImages();
        if(results==null||results.length==0||results[0]==null)throw new IllegalStateException("Java source returned no classification");
        byte[] labels=new byte[results[0].getDataSize()];
        results[0].exportData(0,labels.length,labels);
        Files.write(Paths.get(args[9]),labels);
        System.out.println("completed="+algorithm.isCompleted()+" dimensions="+
                Arrays.toString(results[0].getExtents())+" result_sha256="+sha256(labels));
        image.disposeLocal();
        for(ModelImage result:results)if(result!=null)result.disposeLocal();
    }
}
