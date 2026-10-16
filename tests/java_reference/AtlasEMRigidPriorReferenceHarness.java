import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import java.security.MessageDigest;

/**
 * Small selected-source excerpt for DemonToadDeformableAtlas rigid prior
 * bounding boxes, transformed shape values, and AtlasEM masked normalization.
 * It deliberately does not model registration or descriptor parsing.
 */
public final class AtlasEMRigidPriorReferenceHarness {
    private static final int IX=15, IY=13, IZ=11;
    private static final int TX=13, TY=11, TZ=9, CLASSES=2;
    private static final float[] IMAGE_RES={0.9f,1.1f,1.3f};
    private static final float[] TEMPLATE_RES={1.0f,1.0f,1.0f};
    private static final float[] TRANSFORM={0.07f,-0.04f,0.03f,0.25f,-0.30f,0.15f};

    private static float[][][][] shapes() {
        float[][][][] values=new float[CLASSES][TX][TY][TZ];
        for(int k=0;k<CLASSES;k++) for(int x=0;x<TX;x++) for(int y=0;y<TY;y++) for(int z=0;z<TZ;z++) {
            boolean inside=(x>=2 && x<=10 && y>=1 && y<=9 && z>=1 && z<=7);
            if(inside) {
                int code=(x*37+y*23+z*17+k*29)%101;
                values[k][x][y][z]=0.08f+(float)code/125.0f;
            }
        }
        return values;
    }

    private static float interpolate(float[][][] shape,float x,float y,float z) {
        return ImageFunctions.linearClosestInterpolation(shape,x,y,z,TX,TY,TZ);
    }

    private static int[] bounds(float[][][] shape,float[][] matrix) {
        int[] lo={IX,IY,IZ}, hi={0,0,0};
        for(int x=0;x<IX;x++) for(int y=0;y<IY;y++) for(int z=0;z<IZ;z++) {
            float tx=matrix[0][0]*x+matrix[0][1]*y+matrix[0][2]*z+matrix[0][3];
            float ty=matrix[1][0]*x+matrix[1][1]*y+matrix[1][2]*z+matrix[1][3];
            float tz=matrix[2][0]*x+matrix[2][1]*y+matrix[2][2]*z+matrix[2][3];
            if(interpolate(shape,tx,ty,tz)>0.0f) {
                lo[0]=Math.min(lo[0],x);lo[1]=Math.min(lo[1],y);lo[2]=Math.min(lo[2],z);
                hi[0]=Math.max(hi[0],x);hi[1]=Math.max(hi[1],y);hi[2]=Math.max(hi[2],z);
            }
        }
        return new int[]{lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]};
    }

    public static void main(String[] args) throws Exception {
        ParametricTransform model=new ParametricTransform("rigid",IX/2.0f,IY/2.0f,IZ/2.0f,
                IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],IX,IY,IZ,
                TX/2.0f,TY/2.0f,TZ/2.0f,TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],TX,TY,TZ);
        float[][] matrix=new float[3][4];
        model.precomputeImageToTemplateMatrix(matrix,TRANSFORM,ParametricTransform.computeRotation(TRANSFORM),1.0f);
        float[][][][] shapes=shapes();
        int[][] box=new int[CLASSES][];
        for(int k=0;k<CLASSES;k++) box[k]=bounds(shapes[k],matrix);
        float[][][][] output=new float[IX][IY][IZ][CLASSES];
        int active=0;
        for(int x=0;x<IX;x++) for(int y=0;y<IY;y++) for(int z=0;z<IZ;z++) {
            boolean mask=(x*3+y*5+z*7)%9!=0;
            if(mask) active++;
            for(int k=0;k<CLASSES;k++) {
                int[] b=box[k];
                float value;
                if(x<b[0]||x>=b[3]||y<b[1]||y>=b[4]||z<b[2]||z>=b[5]) value=(k==0)?1.0f:0.0f;
                else {
                    float tx=matrix[0][0]*x+matrix[0][1]*y+matrix[0][2]*z+matrix[0][3];
                    float ty=matrix[1][0]*x+matrix[1][1]*y+matrix[1][2]*z+matrix[1][3];
                    float tz=matrix[2][0]*x+matrix[2][1]*y+matrix[2][2]*z+matrix[2][3];
                    value=Math.max(0.0f,Math.min(1.0f,2.0f*(interpolate(shapes[k],tx,ty,tz)-0.5f)+0.5f));
                }
                output[x][y][z][k]=value;
            }
            if(mask) {
                float norm=output[x][y][z][0]+output[x][y][z][1];
                if(norm!=0.0f) for(int k=0;k<CLASSES;k++) output[x][y][z][k]/=norm;
            }
        }
        MessageDigest digest=MessageDigest.getInstance("SHA-256");
        for(int x=0;x<IX;x++) for(int y=0;y<IY;y++) for(int z=0;z<IZ;z++) for(int k=0;k<CLASSES;k++) {
            int bits=Float.floatToIntBits(output[x][y][z][k]);
            digest.update((byte)(bits>>>24));digest.update((byte)(bits>>>16));digest.update((byte)(bits>>>8));digest.update((byte)bits);
        }
        StringBuilder hex=new StringBuilder();for(byte value:digest.digest()) hex.append(String.format("%02x",value));
        System.out.printf("{\"source\":\"DemonToadDeformableAtlas rigid transformed shape excerpt + AtlasEM computeAtlas normalization\",\"atlas_source_sha256\":\"96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a\",\"parametric_transform_sha256\":\"22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac\",\"image_functions_sha256\":\"be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892\",\"output_sha256\":\"%s\",\"active_voxels\":%d,\"bounds\":[",hex,active);
        for(int k=0;k<CLASSES;k++) { if(k>0)System.out.print(",");System.out.printf("[%d,%d,%d,%d,%d,%d]",box[k][0],box[k][1],box[k][2],box[k][3],box[k][4],box[k][5]); }
        System.out.println("]}");
    }
}
