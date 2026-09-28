/*
 * Test-only formula excerpt for RotationMatrix.derivatives and the RIGID
 * branch of ParametricTransform.imageToTemplateDerivatives in the selected
 * VandyCRUISE source. It uses a fixed nonsingular synthetic rotation/grid.
 */
public final class ToadsRigidJacobianReferenceHarness {
    private static float[][] rotationDerivative(float a, float b, float c,
            float da, float db, float dc) {
        float a2=a*a, b2=b*b, c2=c*c, ab=a*b, ac=a*c, bc=b*c;
        float x2=1.0f-a2-b2-c2, x=(float)Math.sqrt(x2);
        float ax=x*a, bx=x*b, cx=x*c;
        float[][] r=new float[3][3];
        r[0][0]=-4.0f*(b*db+c*dc);
        r[1][0]=2.0f/x*((bx-ac)*da+(ax-bc)*db+(x2-c2)*dc);
        r[2][0]=2.0f/x*((cx+ab)*da+(b2-x2)*db+(ax+bc)*dc);
        r[0][1]=2.0f/x*((bx+ac)*da+(ax+bc)*db+(c2-x2)*dc);
        r[1][1]=-4.0f*(a*da+c*dc);
        r[2][1]=2.0f/x*((x2-a2)*da+(cx-ab)*db+(bx-ac)*dc);
        r[0][2]=2.0f/x*((cx-ab)*da+(x2-b2)*db+(ax-bc)*dc);
        r[1][2]=2.0f/x*((a2-x2)*da+(cx+ab)*db+(bx+ac)*dc);
        r[2][2]=-4.0f*(a*da+b*db);
        return r;
    }

    public static void main(String[] args) {
        float[] t={0.14f,0.07f,-0.05f,1.0f,2.0f,3.0f};
        float[] voxel={2.0f,1.0f,3.0f};
        float[] center={2.5f,2.0f,1.5f};
        float[] imageRes={0.8f,1.2f,2.0f};
        float[] templateRes={0.9f,1.1f,1.3f};
        float[] phys=new float[3];
        for (int i=0;i<3;i++) phys[i]=(voxel[i]-center[i])*imageRes[i];
        float[][][] dR={rotationDerivative(t[0],t[1],t[2],1,0,0),
            rotationDerivative(t[0],t[1],t[2],0,1,0),
            rotationDerivative(t[0],t[1],t[2],0,0,1)};
        float[][] jac=new float[3][6];
        for (int axis=0;axis<3;axis++) for (int p=0;p<3;p++) {
            float v=(dR[p][axis][0]*phys[0]+dR[p][axis][1]*phys[1]+dR[p][axis][2]*phys[2])/templateRes[axis];
            jac[axis][p]=v;
        }
        jac[0][3]=1.0f/templateRes[0];
        jac[1][4]=1.0f/templateRes[1];
        jac[2][5]=1.0f/templateRes[2];
        StringBuilder out=new StringBuilder("{\"jacobian\":[");
        for (int i=0;i<3;i++) {
            if(i>0) out.append(',');
            out.append('[');
            for (int j=0;j<6;j++) { if(j>0) out.append(','); out.append(Float.toString(jac[i][j])); }
            out.append(']');
        }
        System.out.println(out.append("]}"));
    }
}
