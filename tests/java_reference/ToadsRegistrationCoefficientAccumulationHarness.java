/*
 * Test-only source-formula excerpt for one voxel in
 * MultipleToadDeformableAtlas.computeExactRegistrationCoefficients. It validates
 * prior-gradient/Jacobian composition and the cost, norm, gradient, diagonal
 * Hessian, norm-gradient, and norm-Hessian contribution order.
 */
public final class ToadsRegistrationCoefficientAccumulationHarness {
    public static void main(String[] args) {
        float[] mem={0.9f,0.25f,0.75f,0.6f};
        float[] prior={1.0f,0.5f,0.0f,0.8f};
        float[][] priorGrad={{0,0,0},{0.1f,-0.2f,0.3f},{-0.4f,0.5f,0.2f},{0.15f,0.05f,-0.1f}};
        float[][] jac={{1,2,3,4,0,0},{0.5f,-1,0.25f,0,5,0},{-0.75f,0.1f,0.3f,0,0,6}};
        boolean[] registered={true,true,false,true};
        float cost=0, norm=0;
        float[] grad=new float[6], hess=new float[6], ngrad=new float[6], nhess=new float[6];
        for (int k=1;k<4;k++) if (registered[k] && prior[k]>0) {
            float p=prior[k], m=mem[k];
            float[] dp=new float[6];
            for (int i=0;i<6;i++) dp[i]=priorGrad[k][0]*jac[0][i]+priorGrad[k][1]*jac[1][i]+priorGrad[k][2]*jac[2][i];
            cost += m*m*p*p;
            norm += p*p;
            for (int i=0;i<6;i++) {
                grad[i] += 2.0f*m*m*p*dp[i];
                hess[i] += 2.0f*m*m*dp[i]*dp[i];
                ngrad[i] += 2.0f*p*dp[i];
                nhess[i] += 2.0f*dp[i]*dp[i];
            }
        }
        StringBuilder out=new StringBuilder("{\"cost\":").append(Float.toString(cost))
            .append(",\"norm\":").append(Float.toString(norm));
        for (String name : new String[]{"gradient","hessian","norm_gradient","norm_hessian"}) {
            float[] values=name.equals("gradient")?grad:name.equals("hessian")?hess:name.equals("norm_gradient")?ngrad:nhess;
            out.append(",\"").append(name).append("\":[");
            for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.toString(values[i]));}
            out.append(']');
        }
        System.out.println(out.append('}'));
    }
}
