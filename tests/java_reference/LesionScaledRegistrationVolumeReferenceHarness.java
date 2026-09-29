import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Formula excerpt for fixed-transform, scaled LesionToadAtlas volume coefficients. */
public final class LesionScaledRegistrationVolumeReferenceHarness {
    private static final int[] IMAGE_DIMS = {13, 11, 9};
    private static final int[] TEMPLATE_DIMS = {20, 18, 16};
    private static final float[] IMAGE_RES = {1.1f, 1.3f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9f, 1.0f, 1.2f};
    private static final float[] TRANSFORM = {0.05f,-0.03f,0.04f,0.1f,-0.1f,0.2f};
    private static final int CLASSES = 4;
    private static final int SUBSAMPLE = 2;

    private static int[] coarse(int[] dims, int scale) {
        return new int[] {dims[0]/scale,dims[1]/scale,dims[2]/scale};
    }

    private static float[][][] shape(int clazz, int[] dims) {
        float[][][] values = new float[dims[0]][dims[1]][dims[2]];
        for (int x=0; x<dims[0]; x++) for (int y=0; y<dims[1]; y++) for (int z=0; z<dims[2]; z++) {
            if (clazz==0) values[x][y][z]=0.9f;
            else {
                float value=0.1f;
                value += 0.01f*clazz;
                value += 0.013f*x;
                value += 0.019f*y;
                value += 0.027f*z;
                value += 0.001f*x*y;
                values[x][y][z]=value;
            }
        }
        return values;
    }

    private static float membership(int x, int y, int z, int clazz) {
        float value=0.23f;
        value += 0.1f*clazz;
        value += 0.017f*x;
        value += 0.011f*y;
        value += 0.007f*z;
        return value;
    }

    private static String vector(float[] values) {
        StringBuilder out=new StringBuilder("[");
        for (int i=0;i<values.length;i++) { if(i>0) out.append(','); out.append(Float.toString(values[i])); }
        return out.append(']').toString();
    }

    private static String run(int scale) {
        int[] input=coarse(IMAGE_DIMS,scale), template=coarse(TEMPLATE_DIMS,scale);
        float[] imageCenter={IMAGE_DIMS[0]/2.0f,IMAGE_DIMS[1]/2.0f,IMAGE_DIMS[2]/2.0f};
        float[] templateCenter={TEMPLATE_DIMS[0]/2.0f,TEMPLATE_DIMS[1]/2.0f,TEMPLATE_DIMS[2]/2.0f};
        ParametricTransform model=new ParametricTransform("rigid",
            imageCenter[0],imageCenter[1],imageCenter[2],
            IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2],
            templateCenter[0],templateCenter[1],templateCenter[2],
            TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],
            TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
        RotationMatrix rotation=model.computeRotationMatrix(TRANSFORM);
        float[][] rot=rotation.getMatrix();
        float[][] dRa=rotation.derivatives(1.0f,0.0f,0.0f);
        float[][] dRb=rotation.derivatives(0.0f,1.0f,0.0f);
        float[][] dRc=rotation.derivatives(0.0f,0.0f,1.0f);
        float[][][][] shapes=new float[CLASSES][][][];
        for (int k=0;k<CLASSES;k++) shapes[k]=shape(k,template);
        float[][][][] mem=new float[input[0]][input[1]][input[2]][CLASSES];
        for (int x=0;x<input[0];x++) for (int y=0;y<input[1];y++) for (int z=0;z<input[2];z++)
            for (int k=0;k<CLASSES;k++) mem[x][y][z][k]=membership(x,y,z,k);

        float cost=0.0f,norm=0.0f;
        float[] gradient=new float[6],hessian=new float[6];
        int voxels=0,terms=0;
        for (int x=0;x<input[0];x+=SUBSAMPLE) for (int y=0;y<input[1];y+=SUBSAMPLE) for (int z=0;z<input[2];z+=SUBSAMPLE) {
            voxels++;
            for (int k=0;k<CLASSES;k++) if(k!=0) {
                float[] mapped=model.imageToTemplate(x,y,z,TRANSFORM,rot,(float)scale);
                float p=ImageFunctions.linearInterpolation(shapes[k],0.0f,mapped[0],mapped[1],mapped[2],template[0],template[1],template[2]);
                if (p>0.0f) {
                    float dx=ImageFunctions.linearInterpolationXderivative(shapes[k],mapped[0],mapped[1],mapped[2],template[0],template[1],template[2]);
                    float dy=ImageFunctions.linearInterpolationYderivative(shapes[k],mapped[0],mapped[1],mapped[2],template[0],template[1],template[2]);
                    float dz=ImageFunctions.linearInterpolationZderivative(shapes[k],mapped[0],mapped[1],mapped[2],template[0],template[1],template[2]);
                    float[][] dXi=model.imageToTemplateDerivatives(x,y,z,TRANSFORM,rot,dRa,dRb,dRc,(float)scale);
                    float m=mem[x][y][z][k];
                    for (int i=0;i<6;i++) {
                        float dp=dx*dXi[0][i]+dy*dXi[1][i]+dz*dXi[2][i];
                        hessian[i] += 2.0f*m*m*dp*dp;
                        gradient[i] += 2.0f*m*m*p*dp;
                    }
                    cost += m*m*p*p;
                    norm += p*p;
                    terms++;
                }
            }
        }
        float[] normalizedGradient=new float[6],normalizedHessian=new float[6];
        for(int i=0;i<6;i++) { normalizedGradient[i]=gradient[i]/norm; normalizedHessian[i]=hessian[i]/norm; }
        return "{\"scale\":"+scale+",\"input_dims\":"+vector(new float[]{input[0],input[1],input[2]})
            +",\"template_dims\":"+vector(new float[]{template[0],template[1],template[2]})
            +",\"cost\":"+Float.toString(cost)+",\"norm\":"+Float.toString(norm)
            +",\"gradient\":"+vector(gradient)+",\"hessian\":"+vector(hessian)
            +",\"energy\":"+Float.toString(cost/norm)+",\"normalized_gradient\":"+vector(normalizedGradient)
            +",\"normalized_hessian\":"+vector(normalizedHessian)+",\"visited_voxels\":"+voxels
            +",\"contributing_class_terms\":"+terms+'}';
    }

    public static void main(String[] args) {
        System.out.println("{\"full_image_dims\":"+vector(new float[]{IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2]})
            +",\"full_template_dims\":"+vector(new float[]{TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]})
            +",\"image_resolutions\":"+vector(IMAGE_RES)+",\"template_resolutions\":"+vector(TEMPLATE_RES)
            +",\"image_center\":"+vector(new float[]{IMAGE_DIMS[0]/2.0f,IMAGE_DIMS[1]/2.0f,IMAGE_DIMS[2]/2.0f})
            +",\"template_center\":"+vector(new float[]{TEMPLATE_DIMS[0]/2.0f,TEMPLATE_DIMS[1]/2.0f,TEMPLATE_DIMS[2]/2.0f})
            +",\"transform\":"+vector(TRANSFORM)+",\"subsample\":"+SUBSAMPLE
            +",\"classes\":"+CLASSES+",\"cases\":["+run(2)+","+run(4)+"]}");
    }
}
