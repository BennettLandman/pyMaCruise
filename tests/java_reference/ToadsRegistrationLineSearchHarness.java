/*
 * Test-only excerpt of MultipleToadDeformableAtlas.registerGradientDescent.
 * The constant/threshold energies make accept and repeated-rejection
 * behavior deterministic while retaining source float32 update order.
 */
public final class ToadsRegistrationLineSearchHarness {
    private static float[] run(float initialEnergy) {
        float[] transform={0.0f}, gradient={1.0f}, hessian={1.0f};
        float lambda=1.0f, factor=1.5f, energy=initialEnergy;
        int rejected=0;
        boolean accepted=false;
        while (true) {
            float[] trial=new float[1];
            trial[0]=transform[0]+lambda/hessian[0]*gradient[0];
            energy=1.0f;
            if (energy>initialEnergy) {
                transform=trial;
                lambda=lambda*factor;
                accepted=true;
                break;
            } else {
                lambda=lambda/factor;
                rejected++;
                if (rejected>10) break;
            }
        }
        float change=(energy-initialEnergy)/energy;
        return new float[]{transform[0],energy,lambda,change,accepted?1.0f:0.0f,rejected};
    }

    private static void append(String name, float[] values, StringBuilder out) {
        if(out.length()>1) out.append(',');
        out.append('"').append(name).append("\":[");
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.toString(values[i]));}
        out.append(']');
    }

    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{");
        append("accepted",run(0.0f),out);
        append("rejected",run(2.0f),out);
        System.out.println(out.append('}'));
    }
}
