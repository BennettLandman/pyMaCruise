import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Formula excerpt for one coarse voxel of LesionToadAtlas registration. */
public final class LesionScaledRegistrationVoxelReferenceHarness {
    private static final int[] IMAGE_DIMS = {21, 19, 17};
    private static final int[] TEMPLATE_DIMS = {24, 24, 24};
    private static final float[] IMAGE_RES = {1.1f, 1.3f, 1.5f};
    private static final float[] TEMPLATE_RES = {0.9f, 1.0f, 1.2f};
    private static final float[] BASE_TRANSFORM = {0.05f, -0.03f, 0.04f, 0.1f, 0.2f, -0.1f};
    private static final float[] MEMBERSHIPS = {0.9f, 0.5f, 0.75f, 0.6f};

    private static int[] divided(int[] values, int divisor) {
        return new int[] {values[0]/divisor, values[1]/divisor, values[2]/divisor};
    }

    private static float[][][] shape(int clazz, int[] dims) {
        float[][][] out = new float[dims[0]][dims[1]][dims[2]];
        for (int x=0; x<dims[0]; x++) for (int y=0; y<dims[1]; y++) for (int z=0; z<dims[2]; z++) {
            if (clazz==0) out[x][y][z]=0.9f;
            else if (clazz==1) out[x][y][z]=0.1f+0.01f*x+0.02f*y+0.03f*z+0.001f*x*y;
            else if (clazz==2) out[x][y][z]=0.15f+0.008f*x+0.013f*y+0.021f*z+0.0005f*x*z;
            else out[x][y][z]=0.0f;
        }
        return out;
    }

    private static String vec(float[] values) {
        StringBuilder out=new StringBuilder("[");
        for (int i=0;i<values.length;i++) { if(i>0) out.append(','); out.append(Float.toString(values[i])); }
        return out.append(']').toString();
    }

    private static String matrix(float[][] values) {
        StringBuilder out=new StringBuilder("[");
        for (int i=0;i<values.length;i++) { if(i>0) out.append(','); out.append(vec(values[i])); }
        return out.append(']').toString();
    }

    private static void appendCase(StringBuilder output, int scale, boolean outside) {
        int[] imageCoarse=divided(IMAGE_DIMS,scale);
        int[] templateCoarse=divided(TEMPLATE_DIMS,scale);
        int[] point={imageCoarse[0]/2,imageCoarse[1]/2,imageCoarse[2]/2};
        float[] imageCenter={IMAGE_DIMS[0]/2.0f,IMAGE_DIMS[1]/2.0f,IMAGE_DIMS[2]/2.0f};
        float[] templateCenter={8.0f,8.0f,8.0f};
        float[] trans=BASE_TRANSFORM.clone();
        if (outside) trans[3]=-100.0f;
        ParametricTransform model=new ParametricTransform("rigid",
                imageCenter[0],imageCenter[1],imageCenter[2],
                IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2],
                templateCenter[0],templateCenter[1],templateCenter[2],
                TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],
                TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
        RotationMatrix rotation=model.computeRotationMatrix(trans);
        float[][] rot=rotation.getMatrix();
        float[][] dRa=rotation.derivatives(1.0f,0.0f,0.0f);
        float[][] dRb=rotation.derivatives(0.0f,1.0f,0.0f);
        float[][] dRc=rotation.derivatives(0.0f,0.0f,1.0f);
        float[] mapped=model.imageToTemplate(point[0],point[1],point[2],trans,rot,(float)scale);
        float[][] dXi=model.imageToTemplateDerivatives(point[0],point[1],point[2],trans,
                rot,dRa,dRb,dRc,(float)scale);
        float[][] jacobian=new float[3][6];
        for (int axis=0;axis<3;axis++) for (int parameter=0;parameter<6;parameter++)
            jacobian[axis][parameter]=dXi[axis][parameter];

        float[] priors=new float[4];
        float[][] priorGradients=new float[4][3];
        float cost=0.0f,norm=0.0f;
        float[] gradient=new float[6],hessian=new float[6];
        for (int k=0;k<4;k++) {
            float[][][] map=shape(k,templateCoarse);
            priors[k]=k==0?0.0f:ImageFunctions.linearInterpolation(map,0.0f,
                    mapped[0],mapped[1],mapped[2],templateCoarse[0],templateCoarse[1],templateCoarse[2]);
            if (k>0 && priors[k]>0.0f) {
                priorGradients[k][0]=ImageFunctions.linearInterpolationXderivative(map,mapped[0],mapped[1],mapped[2],
                        templateCoarse[0],templateCoarse[1],templateCoarse[2]);
                priorGradients[k][1]=ImageFunctions.linearInterpolationYderivative(map,mapped[0],mapped[1],mapped[2],
                        templateCoarse[0],templateCoarse[1],templateCoarse[2]);
                priorGradients[k][2]=ImageFunctions.linearInterpolationZderivative(map,mapped[0],mapped[1],mapped[2],
                        templateCoarse[0],templateCoarse[1],templateCoarse[2]);
                float[] dp=new float[6];
                for (int i=0;i<6;i++)
                    dp[i]=priorGradients[k][0]*dXi[0][i]+priorGradients[k][1]*dXi[1][i]+priorGradients[k][2]*dXi[2][i];
                float m=MEMBERSHIPS[k],p=priors[k];
                cost += m*m*p*p;
                norm += p*p;
                for (int i=0;i<6;i++) {
                    gradient[i] += 2.0f*m*m*p*dp[i];
                    hessian[i] += 2.0f*m*m*dp[i]*dp[i];
                }
            }
        }
        if (output.charAt(output.length()-1)!='[') output.append(',');
        output.append("{\"scale\":").append(scale)
            .append(",\"outside\":").append(outside)
            .append(",\"point\":").append(vec(new float[]{point[0],point[1],point[2]}))
            .append(",\"image_dims\":").append(vec(new float[]{imageCoarse[0],imageCoarse[1],imageCoarse[2]}))
            .append(",\"template_dims\":").append(vec(new float[]{templateCoarse[0],templateCoarse[1],templateCoarse[2]}))
            .append(",\"image_center\":").append(vec(imageCenter))
            .append(",\"template_center\":").append(vec(templateCenter))
            .append(",\"image_resolutions\":").append(vec(IMAGE_RES))
            .append(",\"template_resolutions\":").append(vec(TEMPLATE_RES))
            .append(",\"transform\":").append(vec(trans))
            .append(",\"mapped\":").append(vec(mapped))
            .append(",\"priors\":").append(vec(priors))
            .append(",\"prior_gradients_xyz\":").append(matrix(priorGradients))
            .append(",\"jacobian\":").append(matrix(jacobian))
            .append(",\"cost\":").append(Float.toString(cost))
            .append(",\"norm\":").append(Float.toString(norm))
            .append(",\"gradient\":").append(vec(gradient))
            .append(",\"hessian\":").append(vec(hessian)).append('}');
    }

    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{\"cases\":[");
        for (int scale : new int[] {2,4,8}) {
            appendCase(out,scale,false);
            appendCase(out,scale,true);
        }
        System.out.println(out.append("]}"));
    }
}
