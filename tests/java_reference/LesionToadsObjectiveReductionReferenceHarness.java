import java.io.File;
import java.util.Locale;
import java.util.Scanner;

/** Selected VandyCRUISE objective accumulator excerpt, driven by NIfTI floats. */
public final class LesionToadsObjectiveReductionReferenceHarness {
    private static final int NX=5, NY=6, NZ=7, NC=3, CHANNELS=11;
    private static long bits(float value){return Float.floatToIntBits(value)&0xffffffffL;}

    public static void main(String[] args)throws Exception{
        Locale.setDefault(Locale.ROOT);
        if(args.length!=1)throw new IllegalArgumentException("pass flat x/y/z/channel float32 text export");
        float[][][][] input=new float[NX][NY][NZ][CHANNELS];
        try(Scanner scan=new Scanner(new File(args[0]))){
            scan.useLocale(Locale.ROOT);
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)input[x][y][z][c]=scan.nextFloat();
        }
        float total=0.0f;
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,6,7],\"classes\":3,\"voxels\":[");
        boolean first=true;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float denominator=input[x][y][z][6];
            if(denominator>0.0f){
                for(int k=0;k<NC;k++){
                    float term=input[x][y][z][3+k]*input[x][y][z][k];
                    term=term*input[x][y][z][k];
                    total=total+term;
                }
            }
            float lesionTerm=input[x][y][z][9]*input[x][y][z][7];
            lesionTerm=lesionTerm*input[x][y][z][7];
            float blackHoleTerm=input[x][y][z][10]*input[x][y][z][8];
            blackHoleTerm=blackHoleTerm*input[x][y][z][8];
            total=total+(lesionTerm+blackHoleTerm);
            if(!first)out.append(',');first=false;
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z)
                .append("],\"running_total_bits\":").append(bits(total)).append('}');
        }
        System.out.println(out.append("],\"total_bits\":").append(bits(total)).append('}').toString());
    }
}
