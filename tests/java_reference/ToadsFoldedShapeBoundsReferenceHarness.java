import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/** Larger folded synthetic support scan using selected VandyCRUISE source APIs. */
public final class ToadsFoldedShapeBoundsReferenceHarness {
    private static final int[] IMAGE_DIMS = {43, 39, 35};
    private static final int[] TEMPLATE_DIMS = {39, 35, 31};
    private static final float[] RES = {1.0f, 1.0f, 1.0f};

    private static final class Case {
        final String name;
        final float[] transform;
        Case(String name, float[] transform) { this.name=name; this.transform=transform; }
    }

    private static float[][][] foldedPrior() {
        float[][][] values = new float[TEMPLATE_DIMS[0]][TEMPLATE_DIMS[1]][TEMPLATE_DIMS[2]];
        double cx=(TEMPLATE_DIMS[0]-1)/2.0, cy=(TEMPLATE_DIMS[1]-1)/2.0;
        double cz=(TEMPLATE_DIMS[2]-1)/2.0;
        for (int x=0;x<TEMPLATE_DIMS[0];x++) for (int y=0;y<TEMPLATE_DIMS[1];y++)
            for (int z=0;z<TEMPLATE_DIMS[2];z++) {
                double dx=x-cx, dy=y-cy, dz=z-cz;
                double radius=Math.sqrt(dx*dx+dy*dy+0.35*dz*dz);
                double fold=2.4*Math.sin(dx*0.48+dz*0.19)
                        +1.7*Math.sin(dy*0.53-dz*0.31)
                        +1.1*Math.sin((dx+dy)*0.29+dz*0.37);
                double distance=Math.abs(radius-(12.0+fold));
                values[x][y][z]=(float)Math.max(0.0, 1.0-distance/2.6);
            }
        return values;
    }

    private static String triple(int[] values) {
        return String.format("[%d,%d,%d]",values[0],values[1],values[2]);
    }

    private static String vector(float[] values) {
        return String.format("[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g]",values[0],values[1],
                values[2],values[3],values[4],values[5]);
    }

    public static void main(String[] args) {
        Case[] cases={
            new Case("identity",new float[]{0,0,0,0,0,0}),
            new Case("translated",new float[]{0,0,0,2.25f,-1.5f,1.0f}),
            new Case("rotated_translated",new float[]{0.11f,-0.07f,0.05f,-1.25f,1.75f,-0.5f})
        };
        float[][][] prior=foldedPrior();
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadDeformableAtlas.computeTransformedShapeBoundingBox excerpt\",");
        System.out.print("\"source_sha256\":\"e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2\",");
        System.out.print("\"parametric_transform_sha256\":\"22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac\",");
        System.out.print("\"image_dims\":[43,39,35],\"template_dims\":[39,35,31],\"resolutions\":[1,1,1],\"cases\":[");
        for(int c=0;c<cases.length;c++) {
            Case test=cases[c];
            float[] centersImage={IMAGE_DIMS[0]/2.0f,IMAGE_DIMS[1]/2.0f,IMAGE_DIMS[2]/2.0f};
            float[] centersTemplate={TEMPLATE_DIMS[0]/2.0f,TEMPLATE_DIMS[1]/2.0f,TEMPLATE_DIMS[2]/2.0f};
            ParametricTransform model=new ParametricTransform("rigid",centersImage[0],centersImage[1],centersImage[2],
                    RES[0],RES[1],RES[2],IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2],
                    centersTemplate[0],centersTemplate[1],centersTemplate[2],
                    RES[0],RES[1],RES[2],TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
            float[][] matrix=new float[3][4];
            float[][] rotation=ParametricTransform.computeRotation(test.transform);
            model.precomputeImageToTemplateMatrix(matrix,test.transform,rotation,1.0f);
            int[] low={IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2]}, high={0,0,0};
            int count=0;
            for(int x=0;x<IMAGE_DIMS[0];x++) for(int y=0;y<IMAGE_DIMS[1];y++) for(int z=0;z<IMAGE_DIMS[2];z++) {
                float tx=matrix[0][0]*x+matrix[0][1]*y+matrix[0][2]*z+matrix[0][3];
                float ty=matrix[1][0]*x+matrix[1][1]*y+matrix[1][2]*z+matrix[1][3];
                float tz=matrix[2][0]*x+matrix[2][1]*y+matrix[2][2]*z+matrix[2][3];
                if(ImageFunctions.linearClosestInterpolation(prior,tx,ty,tz,
                        TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2])>0.0f) {
                    low[0]=Math.min(low[0],x); low[1]=Math.min(low[1],y); low[2]=Math.min(low[2],z);
                    high[0]=Math.max(high[0],x); high[1]=Math.max(high[1],y); high[2]=Math.max(high[2],z);
                    count++;
                }
            }
            if(c>0) System.out.print(',');
            System.out.printf("{\"name\":\"%s\",\"transform\":%s,\"lower\":%s,\"upper\":%s,\"support_count\":%d}",
                    test.name,vector(test.transform),triple(low),triple(high),count);
        }
        System.out.println("]}");
    }
}
