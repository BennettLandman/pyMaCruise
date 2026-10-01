import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import java.util.Locale;

/**
 * Selected-Java reference for the post-pyramid initial-alignment bounds scan.
 *
 * The supplied transform is the final output of the separately tested P39
 * initial-membership/pyramid composition. This harness validates the next
 * wrapper operation, LesionToadAtlas.computeTransformedShapeBoundingBox(),
 * using the selected ParametricTransform and ImageFunctions methods.
 */
public final class LesionInitialAlignmentBoundsReferenceHarness {
    private static final int NX=33, NY=33, NZ=33, K=3;
    private static final int SX=35, SY=37, SZ=36;
    private static final String ALGORITHM_SHA=
        "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final String ATLAS_SHA=
        "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e";
    private static final String TRANSFORM_SHA=
        "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac";
    private static final String ROTATION_SHA=
        "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0";
    private static final String IMAGE_FUNCTIONS_SHA=
        "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892";

    private static float[][][] shape(int k) {
        float[][][] image=new float[SX][SY][SZ];
        for(int x=0;x<SX;x++) for(int y=0;y<SY;y++) for(int z=0;z<SZ;z++) {
            if(x<3+k || x>=31-k || y<4 || y>=33-k || z<2+k || z>=32) continue;
            int fold=(3*x+5*y+7*z+11*x*y+13*y*z+k)%11-5;
            float value=(float)(0.72-0.12*k);
            value=(float)(value+(float)(0.0015*x));
            value=(float)(value+(float)(0.0011*y));
            value=(float)(value+(float)(0.0009*z));
            image[x][y][z]=(float)(value+(float)(fold*0.0002));
        }
        return image;
    }

    private static String hex(float value) {
        return String.format(Locale.ROOT,"%08x",Float.floatToRawIntBits(value));
    }

    private static String triple(int[] values) {
        return "["+values[0]+","+values[1]+","+values[2]+"]";
    }

    public static void main(String[] args) {
        if(args.length!=6) throw new IllegalArgumentException("expected six transform float bit patterns");
        float[] transform=new float[6];
        for(int i=0;i<6;i++) transform[i]=Float.intBitsToFloat((int)Long.parseLong(args[i],16));
        float[] imageCenter={NX/2.0f,NY/2.0f,NZ/2.0f};
        ParametricTransform model=new ParametricTransform("rigid",
            imageCenter[0],imageCenter[1],imageCenter[2],1.0f,1.0f,1.0f,NX,NY,NZ,
            SX/2.0f,SY/2.0f,SZ/2.0f,1.0f,1.0f,1.0f,SX,SY,SZ);
        RotationMatrix rotation=model.computeRotationMatrix(transform);
        float[][] rot=rotation.getMatrix();
        float[][] affine=new float[3][4];
        model.precomputeImageToTemplateMatrix(affine,transform,rot,1.0f);
        StringBuilder out=new StringBuilder("{\"algorithm_source_sha256\":\"").append(ALGORITHM_SHA)
            .append("\",\"atlas_source_sha256\":\"").append(ATLAS_SHA)
            .append("\",\"transform_source_sha256\":\"").append(TRANSFORM_SHA)
            .append("\",\"rotation_source_sha256\":\"").append(ROTATION_SHA)
            .append("\",\"image_functions_source_sha256\":\"").append(IMAGE_FUNCTIONS_SHA)
            .append("\",\"template_support_rule\":\"x=[3+k,31-k),y=[4,33-k),z=[2+k,32)\"")
            .append(",\"transform_bits\":[");
        for(int i=0;i<6;i++){if(i>0)out.append(',');out.append('"').append(hex(transform[i])).append('"');}
        out.append("],\"bounds\":[");
        for(int k=0;k<K;k++) {
            if(k>0)out.append(',');
            float[][][] values=shape(k);
            int[] min={NX,NY,NZ},max={0,0,0};int support=0;
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                float tx=affine[0][0]*x+affine[0][1]*y+affine[0][2]*z+affine[0][3];
                float ty=affine[1][0]*x+affine[1][1]*y+affine[1][2]*z+affine[1][3];
                float tz=affine[2][0]*x+affine[2][1]*y+affine[2][2]*z+affine[2][3];
                if(ImageFunctions.linearInterpolation(values,0.0f,tx,ty,tz,SX,SY,SZ)>0.0f) {
                    support++;
                    min[0]=Math.min(min[0],x);min[1]=Math.min(min[1],y);min[2]=Math.min(min[2],z);
                    max[0]=Math.max(max[0],x);max[1]=Math.max(max[1],y);max[2]=Math.max(max[2],z);
                }
            }
            out.append("{\"class\":").append(k).append(",\"lower_xyz\":").append(triple(min))
                .append(",\"upper_xyz\":").append(triple(max)).append(",\"support_count\":").append(support).append('}');
        }
        out.append("],\"call_order\":[\"atlas.registerShapesPyramid\",\"atlas.computeTransformedShapeBoundingBox\"]}\n");
        System.out.print(out);
    }
}
