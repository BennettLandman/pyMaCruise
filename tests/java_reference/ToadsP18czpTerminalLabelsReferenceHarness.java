import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;

/** Selected updateClassificationFromLabels loop on explicit fixed-stop arrays. */
public final class ToadsP18czpTerminalLabelsReferenceHarness {
    public static void main(String[] args) throws Exception {
        if(args.length!=5)throw new IllegalArgumentException("usage: harness state.bin nx ny nz classes");
        int nx=Integer.parseInt(args[1]),ny=Integer.parseInt(args[2]),nz=Integer.parseInt(args[3]),classes=Integer.parseInt(args[4]);
        int count=nx*ny*nz;
        byte[] payload=Files.readAllBytes(Path.of(args[0]));
        if(payload.length!=count+classes+count)throw new IllegalArgumentException("bad explicit classification payload length");
        ByteBuffer b=ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        byte[] segmentation=new byte[count],template=new byte[classes],previous=new byte[count],output=new byte[count];
        b.get(segmentation);b.get(template);b.get(previous);
        int changed=0,index=0;
        // This is the selected x/y/z order and Java byte assignment/count path.
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){
            int label=segmentation[index];output[index]=template[label];if(previous[index]!=output[index])changed++;index++;
        }
        StringBuilder values=new StringBuilder("[");
        for(int i=0;i<count;i++){if(i>0)values.append(',');values.append(output[i]);}
        System.out.println("{\"reference_scope\":\"selected Java updateClassificationFromLabels excerpt with fixed stop decision; explicit arrays\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"shape\":["+nx+","+ny+","+nz+"],\"classification\":"+values.append(']')+",\"changed_voxels\":"+changed+"}");
    }
}
