/*
 * Test-only deterministic excerpt of the outer x/y/z/class loops in
 * MultipleToadDeformableAtlas.computeExactRegistrationCoefficients. Synthetic
 * arrays supply sampled priors, prior derivatives, memberships, and rigid
 * Jacobians. No MIPAV/JIST image runtime is involved.
 */
public final class ToadsRegistrationFullCoefficientsHarness {
    private static final int NX=3, NY=2, NZ=2, NC=4;
    private static float mem(int x,int y,int z,int k){return ((7*x+3*y+5*z+2*k)%17)/16.0f;}
    private static float prior(int x,int y,int z,int k){return ((2*x+5*y+3*z+4*k)%11)/10.0f;}
    private static float priorGrad(int x,int y,int z,int k,int axis){
        int value=axis==0?x+2*y+3*z+k:axis==1?2*x+y+z+2*k:3*x+y+2*z+k;
        return (value%7-3)/10.0f;
    }
    private static float jac(int x,int y,int z,int axis,int parameter){
        return ((x*3+y*5+z*7+axis*2+parameter*4)%13-6)/10.0f;
    }

    private static float[][] evaluate(int offset,int subsample){
        float cost=0,norm=0;
        float[] grad=new float[6],hess=new float[6],ngrad=new float[6],nhess=new float[6];
        boolean[] registered={true,true,false,true};
        for(int x=offset;x<NX;x+=subsample) for(int y=offset;y<NY;y+=subsample)
            for(int z=offset;z<NZ;z+=subsample) for(int k=1;k<NC;k++)
                if(registered[k]) {
                    float p=prior(x,y,z,k), m=mem(x,y,z,k);
                    if(p>0) {
                        float[] dp=new float[6];
                        for(int t=0;t<6;t++)
                            dp[t]=priorGrad(x,y,z,k,0)*jac(x,y,z,0,t)
                                +priorGrad(x,y,z,k,1)*jac(x,y,z,1,t)
                                +priorGrad(x,y,z,k,2)*jac(x,y,z,2,t);
                        cost+=m*m*p*p;
                        norm+=p*p;
                        for(int t=0;t<6;t++) {
                            grad[t]+=2.0f*m*m*p*dp[t];
                            hess[t]+=2.0f*m*m*dp[t]*dp[t];
                            ngrad[t]+=2.0f*p*dp[t];
                            nhess[t]+=2.0f*dp[t]*dp[t];
                        }
                    }
                }
        return new float[][]{{cost,norm},grad,hess,ngrad,nhess};
    }

    private static void append(String name,float[][] values,StringBuilder out){
        if(out.length()>1)out.append(',');
        out.append('"').append(name).append("\":{");
        String[] keys={"cost_norm","gradient","hessian","norm_gradient","norm_hessian"};
        for(int r=0;r<values.length;r++){
            if(r>0)out.append(',');
            out.append('"').append(keys[r]).append("\":[");
            for(int i=0;i<values[r].length;i++){if(i>0)out.append(',');out.append(Float.toString(values[r][i]));}
            out.append(']');
        }
        out.append('}');
    }
    public static void main(String[] args){
        StringBuilder out=new StringBuilder("{");
        append("full_grid",evaluate(0,1),out);
        append("offset_subsample",evaluate(1,2),out);
        System.out.println(out.append('}'));
    }
}
