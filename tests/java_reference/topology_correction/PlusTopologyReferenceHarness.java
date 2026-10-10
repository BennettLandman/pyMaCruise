import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import edu.vanderbilt.masi.plugins.CRUISE.utilities.MASImorph;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
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

    private static void writeRaw(Path path, int nx, int ny, int nz, float[] values) throws Exception {
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(
                new GZIPOutputStream(new FileOutputStream(path.toFile()),1<<18)))) {
            out.writeInt(nx);out.writeInt(ny);out.writeInt(nz);
            for(float value:values) out.writeFloat(value);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("usage: INPUT_BIN OUTPUT_DIR");
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
                destination,source,"scalar_image",1.0e-4f,3,6,18,"paint_mask",paint,
                0.0f,1.0f,true,"object->background");
        algorithm.runAlgorithm();
        if(!algorithm.isCompleted()) throw new IllegalStateException("Algorithm did not complete");
        Path outputDir=Path.of(args[1]); Files.createDirectories(outputDir);
        float[] corrected=destination.getData();
        writeRaw(outputDir.resolve("topology-java.raw.gz"),nx,ny,nz,corrected);

        // Source-derived MAWM step: MASImorph.Minimize3DFilter uses the
        // center plus six face neighbors; its output border remains zero.
        float[][][] atlas=new float[nx][ny][nz];
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++)
            atlas[x][y][z]=mawm[x+nx*y+nx*ny*z];
        MASImorph morph=new MASImorph(new ImageDataFloat(atlas));
        morph.Minimize3DFilter();
        ImageDataInt eroded=(ImageDataInt)morph.GetOutput();
        float[] erodedFlat=new float[volume.length], refined=corrected.clone();
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) {
            int index=x+nx*y+nx*ny*z;
            erodedFlat[index]=eroded.getInt(x,y,z);
            // Literal Refine_Topo_by_MA.m condition and assignment.
            if(refined[index]<200.0f && erodedFlat[index]>0.0f) refined[index]=255.0f;
        }
        writeRaw(outputDir.resolve("topology-mawm-eroded.raw.gz"),nx,ny,nz,erodedFlat);
        writeRaw(outputDir.resolve("topology-mawm-refined.raw.gz"),nx,ny,nz,refined);
        System.out.println("shape="+nx+"x"+ny+"x"+nz+" completed="+algorithm.isCompleted());
    }
}
