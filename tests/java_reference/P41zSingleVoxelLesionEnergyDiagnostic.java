import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.DataInputStream;
import java.io.FileInputStream;

/** One-target-voxel decomposition of LesionToadSegmentation.computeMemberships. */
public final class P41zSingleVoxelLesionEnergyDiagnostic {
    private static final int C=13, K=10, NC=2;
    private static float bits(float x) { return x; }
    private static String word(float x) { return String.format("%08x", Float.floatToIntBits(x)); }
    public static void main(String[] args) throws Exception {
        float[] image=new float[NC], mu=new float[NC*C], sd=new float[NC*C];
        float[] lesionMu=new float[NC], lesionSd=new float[NC], weights=new float[NC];
        float[] shapes=new float[C], shapePriors=new float[C], smoothingPriors=new float[C];
        float[] neighborMemberships=new float[6*C], lesionFaces=new float[6], blackFaces=new float[6];
        float exponent, smoothing, priorCoefficient;
        try (DataInputStream in=new DataInputStream(new FileInputStream(args[0]))) {
            for(int i=0;i<NC;i++) image[i]=in.readFloat();
            for(int i=0;i<NC*C;i++) mu[i]=in.readFloat();
            for(int i=0;i<NC*C;i++) sd[i]=in.readFloat();
            for(int i=0;i<NC;i++) lesionMu[i]=in.readFloat();
            for(int i=0;i<NC;i++) lesionSd[i]=in.readFloat();
            for(int i=0;i<NC;i++) weights[i]=in.readFloat();
            for(int i=0;i<C;i++) shapes[i]=in.readFloat();
            for(int i=0;i<C;i++) shapePriors[i]=in.readFloat();
            for(int i=0;i<C;i++) smoothingPriors[i]=in.readFloat();
            for(int i=0;i<6*C;i++) neighborMemberships[i]=in.readFloat();
            for(int i=0;i<6;i++) lesionFaces[i]=in.readFloat();
            for(int i=0;i<6;i++) blackFaces[i]=in.readFloat();
            exponent=in.readFloat(); smoothing=in.readFloat(); priorCoefficient=in.readFloat();
            if(in.read()!=-1) throw new IllegalArgumentException("trailing bytes");
        }

        // Source raw lesion term: FLAIR residual squared, then PowerTable.lookup.
        float delta=image[1]-lesionMu[1];
        float residualSquared=delta*delta;
        PowerTable table=new PowerTable(0.0f,4.0f,0.000001f,exponent);
        float raw=(float)table.lookup((double)residualSquared,(double)exponent);

        float other=0.0f, lesionSq=0.0f, blackSq=0.0f, self=0.0f;
        for(int f=0;f<6;f++) {
            for(int m=0;m<C;m++) if(m!=K) {
                float q=neighborMemberships[f*C+m]; other+=q*q;
            }
            float l=lesionFaces[f], b=blackFaces[f];
            lesionSq+=l*l; blackSq+=b*b;
        }
        for(int f : new int[]{1,0,3,2,5,4}) {
            float q=neighborMemberships[f*C+K]; self+=q*q;
        }
        float smoothSum=other+self+blackSq;
        float smoothingIncrement=0.5f*smoothingPriors[K]*smoothing*smoothSum;

        float invpriorL=0.0f;
        for(int m=0;m<C;m++) {
            float distance=0.0f;
            for(int c=0;c<NC;c++) {
                float q=lesionMu[c]-mu[c*C+m];
                float weighted=weights[c]*q; weighted=weighted*q;
                float den=sd[c*C+m]+lesionSd[c];
                distance+=(weighted/den)/den;
            }
            float similarity=(m==K)?0.0f:1.0f/(1.0f+((1.0f-0.1f)/0.1f)*distance);
            float p=shapePriors[m]*similarity;
            p=p*shapes[m]; p=p*shapes[m]; invpriorL+=p;
        }
        float priorIncrement=0.5f*priorCoefficient*invpriorL;
        float finalEnergy=raw;
        finalEnergy+=smoothingIncrement;
        finalEnergy+=priorIncrement;
        System.out.printf("{\"shape10\":%.9g,\"residual_squared\":%.9g,\"raw_numL\":%.9g,\"raw_bits\":\"%s\",\"other_neighbor_squares\":%.9g,\"lesion_neighbor_squares\":%.9g,\"black_neighbor_squares\":%.9g,\"class10_neighbor_squares\":%.9g,\"smooth_sum\":%.9g,\"smoothing_increment\":%.9g,\"smoothing_bits\":\"%s\",\"invpriorL\":%.9g,\"prior_increment\":%.9g,\"prior_bits\":\"%s\",\"final_energy\":%.9g,\"final_bits\":\"%s\"}%n",
            shapes[K],residualSquared,raw,word(raw),other,lesionSq,blackSq,self,smoothSum,smoothingIncrement,word(smoothingIncrement),invpriorL,priorIncrement,word(priorIncrement),finalEnergy,word(finalEnergy));
    }
}
