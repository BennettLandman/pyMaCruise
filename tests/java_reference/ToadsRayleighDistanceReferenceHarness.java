import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.GZIPInputStream;

/** Selected MultipleToadSegmentation.rayleighDistance Java source excerpt. */
public final class ToadsRayleighDistanceReferenceHarness {
    private static final String SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final float LOGPI2MPI2=(float)Math.log(Math.PI/2.0f)-(float)Math.PI/2.0f;

    private static float rayleighDistance(float val, float sig) {
        if (val<=0) return 1.0f;
        else return val*val/sig*sig - 2.0f*(float)Math.log(val/sig) + LOGPI2MPI2;
    }

    private static byte[] readGzip(String path) throws Exception {
        ByteArrayOutputStream out=new ByteArrayOutputStream();
        try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){
            byte[] buffer=new byte[4096]; int n;
            while((n=in.read(buffer))>=0) out.write(buffer,0,n);
        }
        return out.toByteArray();
    }

    public static void main(String[] args) throws Exception {
        byte[] raw=readGzip(args[0]);
        ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        int count=b.getShort(42);
        if(b.getInt(0)!=348 || b.getShort(70)!=16 || b.getShort(72)!=32 || b.getShort(40)!=4 ||
                count!=17 || b.getShort(44)!=1 || b.getShort(46)!=1 || b.getShort(48)!=2)
            throw new IllegalArgumentException("unexpected Rayleigh NIfTI fixture header");
        int offset=(int)b.getFloat(108);
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"cases\":[");
        for(int i=0;i<count;i++){
            if(i>0) System.out.print(',');
            // NIfTI stores the last (channel) axis as a contiguous plane.
            float val=b.getFloat(offset+4*i);
            float sig=b.getFloat(offset+4*(count+i));
            float output=rayleighDistance(val,sig);
            System.out.print("{\"val_bits\":"+Float.floatToIntBits(val)+
                    ",\"sig_bits\":"+Float.floatToIntBits(sig)+
                    ",\"output_bits\":"+Float.floatToIntBits(output)+"}");
        }
        System.out.println("]}");
    }
}
