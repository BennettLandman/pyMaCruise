/*
 * Test-only transcription of the RIGID branches of
 * ParametricTransform.imageToTemplate and imageToTemplateDerivatives.
 * Values are evaluated at pyramid scales 1, 2, and 4 with full-resolution
 * centers and anisotropic image/template resolutions.
 */
public final class ToadsScaledRigidReferenceHarness {
    private static float[][] rotation(float a, float b, float c) {
        float norm=(float)Math.sqrt(a*a+b*b+c*c);
        if(norm>1.0f) {
            float magnitude=norm%2.0f;
            if(magnitude>1.0f) {
                a=magnitude*a/norm-2.0f*a/norm;
                b=magnitude*b/norm-2.0f*b/norm;
                c=magnitude*c/norm-2.0f*c/norm;
            } else { a=a*magnitude/norm; b=b*magnitude/norm; c=c*magnitude/norm; }
        }
        float a2=a*a,b2=b*b,c2=c*c,ab=a*b,ac=a*c,bc=b*c;
        float w=(float)Math.sqrt(1.0f-a2-b2-c2);
        return new float[][] {
            {1.0f-2.0f*b2-2.0f*c2,2.0f*ab-2.0f*w*c,2.0f*ac+2.0f*w*b},
            {2.0f*ab+2.0f*w*c,1.0f-2.0f*a2-2.0f*c2,2.0f*bc-2.0f*w*a},
            {2.0f*ac-2.0f*w*b,2.0f*bc+2.0f*w*a,1.0f-2.0f*a2-2.0f*b2}
        };
    }
    private static float[][] derivative(float a,float b,float c,float da,float db,float dc) {
        float a2=a*a,b2=b*b,c2=c*c,ab=a*b,ac=a*c,bc=b*c;
        float x2=1.0f-a2-b2-c2,x=(float)Math.sqrt(x2),ax=x*a,bx=x*b,cx=x*c;
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
    private static void appendArray(StringBuilder out,float[] a) {
        out.append('[');
        for(int i=0;i<a.length;i++) { if(i>0) out.append(','); out.append(Float.toString(a[i])); }
        out.append(']');
    }
    private static void appendMatrix(StringBuilder out,float[][] a) {
        out.append('[');
        for(int i=0;i<a.length;i++) { if(i>0) out.append(','); appendArray(out,a[i]); }
        out.append(']');
    }
    public static void main(String[] args) {
        float[] t={0.14f,0.07f,-0.05f,1.0f,2.0f,3.0f};
        float[][] rot=rotation(t[0],t[1],t[2]);
        float[][][] dr={derivative(t[0],t[1],t[2],1,0,0),derivative(t[0],t[1],t[2],0,1,0),derivative(t[0],t[1],t[2],0,0,1)};
        float[] imageRes={0.8f,1.2f,2.0f},templateRes={0.9f,1.1f,1.3f};
        float[] imageCenter={4.5f,4.0f,3.5f},templateCenter={5.5f,5.0f,4.5f};
        float[] point={1.0f,1.0f,0.0f};
        StringBuilder out=new StringBuilder("{\"cases\":[");
        float[] scales={1.0f,2.0f,4.0f};
        for(int s=0;s<scales.length;s++) {
            if(s>0) out.append(',');
            float scale=scales[s],physical[]=new float[3],mapped[]=new float[3];
            for(int a=0;a<3;a++) physical[a]=(point[a]*scale-imageCenter[a])*imageRes[a];
            for(int a=0;a<3;a++) {
                float v=rot[a][0]*physical[0]+rot[a][1]*physical[1]+rot[a][2]*physical[2]+t[a+3];
                mapped[a]=(v/templateRes[a]+templateCenter[a])/scale;
            }
            float[][] jac=new float[3][6];
            for(int a=0;a<3;a++) for(int p=0;p<3;p++) {
                float v=(dr[p][a][0]*physical[0]+dr[p][a][1]*physical[1]+dr[p][a][2]*physical[2])/templateRes[a];
                jac[a][p]=v/scale;
            }
            jac[0][3]=1.0f/templateRes[0]/scale; jac[1][4]=1.0f/templateRes[1]/scale; jac[2][5]=1.0f/templateRes[2]/scale;
            out.append("{\"scale\":").append(Float.toString(scale)).append(",\"mapped\":"); appendArray(out,mapped);
            out.append(",\"jacobian\":"); appendMatrix(out,jac); out.append('}');
        }
        System.out.println(out.append("]}").toString());
    }
}
