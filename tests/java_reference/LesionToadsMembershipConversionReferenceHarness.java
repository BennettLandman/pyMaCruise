import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected class-energy inversion and direct/prior membership formula. */
public final class LesionToadsMembershipConversionReferenceHarness {
    private static final float ZERO=1e-10f,INF=1e10f;
    private static final float[] ENERGIES={
        Math.nextAfter(ZERO,0.0f),ZERO,Math.nextAfter(ZERO,Float.POSITIVE_INFINITY),
        0.25f,0.25f,0.25f,0.25f,0.25f,0.0f,-0.25f,0.25f
    };
    private static final float[] SHAPES={0.7f,0.7f,0.7f,0.25f,0.2f,0.5f,0.75f,0.0f,0.5f,0.5f,0.5f};
    private static final boolean[] PRIOR={false,false,false,false,true,true,true,true,true,true,false};
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static void run(int i,boolean first) {
        float energy=ENERGIES[i],shape=SHAPES[i];
        float membership=(energy>ZERO)?1.0f/energy:INF;
        if(PRIOR[i])membership*=Numerics.min(1.0f,2.0f*shape);
        if(!first)System.out.print(',');
        System.out.print("{\"energy_bits\":"+bits(energy)+",\"shape_bits\":"+bits(shape)+",\"relation\":\""+(PRIOR[i]?"prior":"direct")+"\",\"membership_bits\":"+bits(membership)+"}");
    }
    public static void main(String[] args) {
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships membership inversion + shapeFactor\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"numerics_sha256\":\"13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb\",\"cases\":[");
        for(int i=0;i<ENERGIES.length;i++)run(i,i==0);
        System.out.println("]}");
    }
}
