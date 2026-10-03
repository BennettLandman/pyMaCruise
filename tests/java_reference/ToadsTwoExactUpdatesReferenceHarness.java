import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;

/**
 * Exercises two source-ordered exact membership updates on the startup
 * phantom. The second update refreshes deviations/similarity, then restores
 * the fixed atlas Intensity Variance Factors before membership evaluation.
 */
public final class ToadsTwoExactUpdatesReferenceHarness {
    private static final int C=11;
    private static final float[][] FACTORS={{.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},{1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}};

    private static Field field(Class<?> type,String name)throws Exception{Field f=type.getDeclaredField(name);f.setAccessible(true);return f;}
    private static Object get(Object target,String name)throws Exception{return field(target.getClass(),name).get(target);}
    private static void set(Object target,String name,Object value)throws Exception{field(target.getClass(),name).set(target,value);}
    private static Object getStatic(String name)throws Exception{Field f=field(ToadsInitialThinningStartupReferenceHarness.class,name);return f.get(null);}
    private static String matrixBits(float[][] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append('[');for(int j=0;j<values[i].length;j++){if(j>0)s.append(',');s.append(Float.floatToIntBits(values[i][j]));}s.append(']');}return s.append(']').toString();}

    public static void main(String[] args)throws Exception{
        if(args.length!=1)throw new IllegalArgumentException("fixture directory required");
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        Constructor<?> ctor=ToadsInitialThinningStartupReferenceHarness.class.getDeclaredConstructor(String.class);ctor.setAccessible(true);
        Object run=ctor.newInstance(args[0]);
        float first=(float)get(run,"objective");
        float[][][][] maps=(float[][][][])get(run,"memberships");
        set(run,"initial",maps);set(run,"objective",0.0f);

        Method refresh=ToadsInitialStatisticsMembershipContinuationReferenceHarness.class.getDeclaredMethod("refreshStatistics",Object.class);refresh.setAccessible(true);
        refresh.invoke(null,run);
        float[][] deviations=(float[][])getStatic("VARIANCE");
        float[][] savedDeviations=new float[2][C];
        for(int c=0;c<2;c++)System.arraycopy(deviations[c],0,savedDeviations[c],0,C);
        for(int c=0;c<2;c++)System.arraycopy(FACTORS[c],0,deviations[c],0,C);
        Method memberships=ToadsInitialThinningStartupReferenceHarness.class.getDeclaredMethod("computeMemberships");memberships.setAccessible(true);memberships.invoke(run);
        float second=(float)get(run,"objective");
        float difference=first-second;
        float numerator=2.0f*difference;
        float denominator=first+second;
        float ratio=Math.abs(numerator/denominator);
        System.out.print("{\"source_order\":[\"first computeExactMemberships\",\"computeExactCentroids\",\"computeExactVariances\",\"computeCentroidSimilarity\",\"second computeExactMemberships\",\"energyRatio\"]");
        System.out.print(",\"first_objective_bits\":"+Float.floatToIntBits(first)+",\"second_objective_bits\":"+Float.floatToIntBits(second)+",\"energy_ratio_bits\":"+Float.floatToIntBits(ratio));
        System.out.print(",\"fixed_intensity_variance_factor_bits\":"+matrixBits(FACTORS)+",\"refreshed_standard_deviation_bits\":"+matrixBits(savedDeviations)+"}");
    }
}
