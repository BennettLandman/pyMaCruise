import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.Locale;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Dynamic excerpt of AlgorithmLesionToads energy-ratio update for two floats. */
public final class LesionDynamicEnergyRatioReferenceHarness {
    private static int bits(float value) { return Float.floatToRawIntBits(value); }
    private static String hex(float value) { return String.format(Locale.ROOT,"%08x",bits(value)); }
    public static void main(String[] args)throws Exception {
        if(args.length!=1)throw new IllegalArgumentException("binary input required");
        float previous,current;
        try(DataInputStream in=new DataInputStream(new FileInputStream(args[0]))) {
            previous=in.readFloat();current=in.readFloat();if(in.read()!=-1)throw new IllegalArgumentException("trailing bytes");
        }
        float difference=previous-current;
        float numerator=2.0f*difference;
        float denominator=previous+current;
        float quotient=numerator/denominator;
        float ratio=Numerics.abs(quotient);
        System.out.print("{\"previous_energy_bits\":\""+hex(previous)+"\",\"current_energy_bits\":\""+hex(current)+"\",\"difference_bits\":\""+hex(difference)+"\",\"numerator_bits\":\""+hex(numerator)+"\",\"denominator_bits\":\""+hex(denominator)+"\",\"quotient_bits\":\""+hex(quotient)+"\",\"energy_ratio_bits\":\""+hex(ratio)+"\"}");
    }
}
