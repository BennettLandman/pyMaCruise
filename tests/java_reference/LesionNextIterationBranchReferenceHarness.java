import edu.jhmi.rad.medic.utilities.Numerics;
import java.io.DataInputStream;
import java.io.FileInputStream;

/** Dynamic excerpt of AlgorithmLesionToads' post-registration control branch. */
public final class LesionNextIterationBranchReferenceHarness {
    public static void main(String[] args) throws Exception {
        if(args.length!=1)throw new IllegalArgumentException("expected binary control inputs");
        float previous,current,threshold;int iteration,nIterMax;boolean narrow;
        try(DataInputStream in=new DataInputStream(new FileInputStream(args[0]))) {
            previous=in.readFloat();current=in.readFloat();threshold=in.readFloat();
            iteration=in.readInt();nIterMax=in.readInt();narrow=in.readBoolean();
            if(in.read()!=-1)throw new IllegalArgumentException("trailing control bytes");
        }
        float ratio=Numerics.abs(2.0f*(previous-current)/(previous+current));
        boolean convergenceStop=ratio<threshold;
        boolean iterationLimitStop=iteration>=Math.abs(nIterMax);
        boolean stop=convergenceStop||iterationLimitStop;
        boolean runGrowing=!narrow&&(!stop||nIterMax>0);
        System.out.print("{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\"");
        System.out.print(",\"ratio_bits\":\""+String.format("%08x",Float.floatToRawIntBits(ratio))+"\"");
        System.out.print(",\"convergence_stop\":"+convergenceStop+",\"iteration_limit_stop\":"+iterationLimitStop);
        System.out.print(",\"stop\":"+stop+",\"run_growing\":"+runGrowing+",\"refresh_relations\":"+(!stop)+"}");
    }
}
