import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** Selected default-atlas generateCruiseOutputs method on supplied synthetic P18 state. */
public final class ToadsP18czpCruiseOutputsReferenceHarness {
    static final String[] NAMES={"Background","SulcalCSF","Ventricles","CerebellarGM","CerebralGM","Caudate","Thalamus","Putamen","Brainstem","CerebellarWM","CerebralWM"};
    static String hash(float[] a)throws Exception{
        MessageDigest d=MessageDigest.getInstance("SHA-256");
        ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        for(float v:a){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}
        StringBuilder s=new StringBuilder();for(byte v:d.digest())s.append(String.format("%02x",v));return s.toString();
    }
    public static void main(String[] args)throws Exception{
        if(args.length!=4)throw new IllegalArgumentException("usage: harness state.bin nx ny nz");
        int nx=Integer.parseInt(args[1]),ny=Integer.parseInt(args[2]),nz=Integer.parseInt(args[3]),classes=NAMES.length,v=nx*ny*nz;
        byte[] bytes=Files.readAllBytes(Path.of(args[0]));
        if(bytes.length!=v+v*classes*4)throw new IllegalArgumentException("bad explicit-output payload length");
        ByteBuffer b=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        byte[] segmentation=new byte[v];b.get(segmentation);float[] memberships=new float[v*classes];for(int i=0;i<memberships.length;i++)memberships[i]=b.getFloat();
        float[][] out=new float[4][v];
        // The selected method's ordered atlas-name dispatch: last matching CSF/GM assignment wins.
        for(int k=0;k<classes;k++)if(NAMES[k].equals("SulcalCSF")||NAMES[k].equals("Sulcal-CSF"))
            for(int p=0;p<v;p++)out[0][p]=memberships[p*classes+k];
        for(int k=0;k<classes;k++)if(NAMES[k].equals("CorticalGM")||NAMES[k].equals("Cerebrum-GM")||NAMES[k].equals("CerebralGM"))
            for(int p=0;p<v;p++)out[1][p]=memberships[p*classes+k];
        for(int k=0;k<classes;k++){
            String name=NAMES[k];
            if(name.equals("SulcalCSF")||name.equals("Sulcal-CSF")||name.equals("CorticalGM")||name.equals("Cerebrum-GM")||name.equals("CerebralGM")||name.startsWith("Cerebell")||name.equals("Brainstem")||name.equals("Background"))continue;
            for(int p=0;p<v;p++){
                out[2][p]=out[2][p]+memberships[p*classes+k];
                if(segmentation[p]==k)out[3][p]=1.0f;
            }
        }
        System.out.println("{\"reference_scope\":\"selected Java MultipleToadSegmentation.generateCruiseOutputs exact-membership path on supplied P18 arrays\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"wrapper_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"shape\":["+nx+","+ny+","+nz+"],\"output_hashes\":[\""+hash(out[0])+"\",\""+hash(out[1])+"\",\""+hash(out[2])+"\",\""+hash(out[3])+"\"]}");
    }
}
