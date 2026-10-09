import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import gov.nih.mipav.model.structures.ModelImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.zip.GZIPOutputStream;

/** Selected AlgorithmTopologyCorrection on held-fixed MaCRUISEplus tissue-stage inputs. */
public final class PlusTopologyReferenceHarness {
    private static float[] readArray(DataInputStream in, int count) throws Exception {
        float[] values = new float[count];
        for (int i=0;i<count;i++) values[i]=in.readFloat();
        return values;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("usage: INPUT_BIN OUTPUT_RAW_GZ");
        int nx,ny,nz; float[] volume,paintValues,mawm;
        try (DataInputStream in=new DataInputStream(new BufferedInputStream(new FileInputStream(args[0]),1<<20))) {
            nx=in.readInt(); ny=in.readInt(); nz=in.readInt(); int count=nx*ny*nz;
            volume=readArray(in,count); paintValues=readArray(in,count); mawm=readArray(in,count);
            if(in.read()!=-1) throw new IllegalArgumentException("unexpected trailing input bytes");
        }
        ModelImage source=new ModelImage(new int[]{nx,ny,nz},"plus_topology_input");
        ModelImage destination=new ModelImage(new int[]{nx,ny,nz},"plus_topology_output");
        source.setData(volume);
        BitSet paint=new BitSet(volume.length);
        for(int i=0;i<paintValues.length;i++) if(paintValues[i]>0.5f) paint.set(i);
        AlgorithmTopologyCorrection algorithm=new AlgorithmTopologyCorrection(
                destination,source,"scalar_image",0.05f,3,6,18,"paint_mask",paint,
                0.0f,1.0f,true,"object->background");
        algorithm.runAlgorithm();
        if(!algorithm.isCompleted()) throw new IllegalStateException("Algorithm did not complete");
        Files.createDirectories(Path.of(args[1]).toAbsolutePath().getParent());
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(
                new GZIPOutputStream(new FileOutputStream(args[1]),1<<18)))) {
            out.writeInt(nx);out.writeInt(ny);out.writeInt(nz);
            for(float value:destination.getData()) out.writeFloat(value);
        }
        System.out.println("shape="+nx+"x"+ny+"x"+nz+" completed="+algorithm.isCompleted());
    }
}
