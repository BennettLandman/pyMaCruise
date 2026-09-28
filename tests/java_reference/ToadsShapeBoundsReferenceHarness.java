import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/** Synthetic support-bound scan using the selected Java transform and sampler APIs. */
public final class ToadsShapeBoundsReferenceHarness {
    private static final int[] IMAGE_DIMS = {6, 5, 4};
    private static final int[] TEMPLATE_DIMS = {6, 5, 4};
    private static final float[] IMAGE_RES = {1.0f, 1.25f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9375f, 0.9375f, 0.9375f};

    private static final class Case {
        final String name;
        final String shapeKind;
        final float[] transform;
        Case(String name, String shapeKind, float[] transform) {
            this.name = name; this.shapeKind = shapeKind; this.transform = transform;
        }
    }

    private static float[][][] shape(String kind) {
        float[][][] values = new float[TEMPLATE_DIMS[0]][TEMPLATE_DIMS[1]][TEMPLATE_DIMS[2]];
        for (int x=0;x<TEMPLATE_DIMS[0];x++) for (int y=0;y<TEMPLATE_DIMS[1];y++)
            for (int z=0;z<TEMPLATE_DIMS[2];z++) {
                if (kind.equals("zero")) values[x][y][z] = 0.0f;
                else if (kind.equals("positive")) values[x][y][z] = 0.25f;
                else values[x][y][z] = x>=1 && x<=4 && y>=1 && y<=3 && z>=1 && z<=2
                        ? (0.1f + 0.03f*x + 0.02f*y + 0.01f*z) : -0.2f;
            }
        return values;
    }

    private static String triple(int[] values) {
        return String.format("[%d,%d,%d]", values[0], values[1], values[2]);
    }

    private static String floatArray(float[] values) {
        return String.format("[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g]", values[0], values[1],
                values[2], values[3], values[4], values[5]);
    }

    public static void main(String[] args) {
        Case[] cases = {
            new Case("identity_box", "box", new float[] {0,0,0,0,0,0}),
            new Case("rigid_anisotropic_box", "box", new float[] {0.12f,-0.08f,0.04f,0.35f,-0.2f,0.1f}),
            new Case("translated_box", "box", new float[] {0,0,0,-1.4f,0.8f,0.4f}),
            new Case("all_positive_touches_all_faces", "positive", new float[] {0,0,0,0,0,0}),
            new Case("no_positive_support", "zero", new float[] {0,0,0,0,0,0})
        };
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadDeformableAtlas.computeTransformedShapeBoundingBox excerpt\",");
        System.out.print("\"source_sha256\":\"e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2\",");
        System.out.print("\"parametric_transform_sha256\":\"22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac\",");
        System.out.print("\"image_dims\":[6,5,4],\"template_dims\":[6,5,4],");
        System.out.print("\"image_resolutions\":[1.0,1.25,1.5],\"template_resolutions\":[0.9375,0.9375,0.9375],\"cases\":[");
        for (int c=0;c<cases.length;c++) {
            Case test=cases[c];
            float[][][] prior=shape(test.shapeKind);
            float x0i=IMAGE_DIMS[0]/2.0f,y0i=IMAGE_DIMS[1]/2.0f,z0i=IMAGE_DIMS[2]/2.0f;
            float x0t=TEMPLATE_DIMS[0]/2.0f,y0t=TEMPLATE_DIMS[1]/2.0f,z0t=TEMPLATE_DIMS[2]/2.0f;
            ParametricTransform model=new ParametricTransform("rigid",x0i,y0i,z0i,
                    IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2],
                    x0t,y0t,z0t,TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],
                    TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
            float[][] matrix=new float[3][4];
            float[][] rotation=ParametricTransform.computeRotation(test.transform);
            model.precomputeImageToTemplateMatrix(matrix,test.transform,rotation,1.0f);
            int[] low={IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2]}, high={0,0,0};
            StringBuilder support=new StringBuilder();
            int count=0;
            for(int x=0;x<IMAGE_DIMS[0];x++) for(int y=0;y<IMAGE_DIMS[1];y++) for(int z=0;z<IMAGE_DIMS[2];z++) {
                float tx=matrix[0][0]*x+matrix[0][1]*y+matrix[0][2]*z+matrix[0][3];
                float ty=matrix[1][0]*x+matrix[1][1]*y+matrix[1][2]*z+matrix[1][3];
                float tz=matrix[2][0]*x+matrix[2][1]*y+matrix[2][2]*z+matrix[2][3];
                float value=ImageFunctions.linearClosestInterpolation(prior,tx,ty,tz,
                        TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
                if(value>0.0f) {
                    low[0]=Math.min(low[0],x); low[1]=Math.min(low[1],y); low[2]=Math.min(low[2],z);
                    high[0]=Math.max(high[0],x); high[1]=Math.max(high[1],y); high[2]=Math.max(high[2],z);
                    if(count++>0) support.append(',');
                    support.append(String.format("[%d,%d,%d]",x,y,z));
                }
            }
            if(c>0) System.out.print(',');
            System.out.printf("{\"name\":\"%s\",\"shape_kind\":\"%s\",\"transform\":%s,\"lower\":%s,\"upper\":%s,\"support_count\":%d,\"support\":[%s]}",
                    test.name,test.shapeKind,floatArray(test.transform),triple(low),triple(high),count,support.toString());
        }
        System.out.println("]}");
    }
}
