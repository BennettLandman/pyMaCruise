/*
 * Standalone formula excerpt for the selected VandyCRUISE rigid pyramid
 * registration. It combines RIGID imageToTemplate/Jacobian math,
 * ImageFunctions closest interpolation and prior derivatives, and the
 * x/y/z/class accumulation in computeExactScaledRegistrationEnergy and
 * computeExactScaledRegistrationCoefficients. No MIPAV/JIST runtime is used.
 */
public final class ToadsScaledRegistrationKernelsReferenceHarness {
    private static final int[] FULL_IMAGE={13,12,11}, FULL_TEMPLATE={17,16,15};
    private static final float[] IMAGE_RES={0.8f,1.2f,1.5f};
    private static final float[] TEMPLATE_RES={0.9f,1.1f,1.3f};
    private static final float[] TRANSFORM={0.14f,0.07f,-0.05f,0.2f,-0.1f,0.3f};
    private static float[][] rotation(float a,float b,float c) {
        double x,a2,b2,c2,an,bn,cn,n,m;
        n=Math.sqrt(a*a+b*b+c*c);
        if(n>1.0) {
            m=n%2.0;
            if(m>1.0) { an=m*a/n-2*a/n; bn=m*b/n-2*b/n; cn=m*c/n-2*c/n; }
            else { an=a*m/n; bn=b*m/n; cn=c*m/n; }
        } else { an=a; bn=b; cn=c; }
        a2=an*an; b2=bn*bn; c2=cn*cn;
        x=1.0-a2-b2-c2; if(x<0.0&&x>-1e-10)x=0.0; x=Math.sqrt(x);
        return new float[][]{{(float)(1.0-2.0*b2-2.0*c2),(float)(2.0*an*bn-2.0*x*cn),(float)(2.0*an*cn+2.0*x*bn)},
            {(float)(2.0*an*bn+2.0*x*cn),(float)(1.0-2.0*a2-2.0*c2),(float)(2.0*bn*cn-2.0*x*an)},
            {(float)(2.0*an*cn-2.0*x*bn),(float)(2.0*bn*cn+2.0*x*an),(float)(1.0-2.0*a2-2.0*b2)}};
    }
    private static float[][] derivative(float a,float b,float c,float da,float db,float dc) {
        float a2=a*a,b2=b*b,c2=c*c,ab=a*b,ac=a*c,bc=b*c;
        float x2=1-a2-b2-c2,x=(float)Math.sqrt(x2),ax=x*a,bx=x*b,cx=x*c;
        float[][] r=new float[3][3];
        r[0][0]=-4*(b*db+c*dc);
        r[1][0]=2/x*((bx-ac)*da+(ax-bc)*db+(x2-c2)*dc);
        r[2][0]=2/x*((cx+ab)*da+(b2-x2)*db+(ax+bc)*dc);
        r[0][1]=2/x*((bx+ac)*da+(ax+bc)*db+(c2-x2)*dc);
        r[1][1]=-4*(a*da+c*dc);
        r[2][1]=2/x*((x2-a2)*da+(cx-ab)*db+(bx-ac)*dc);
        r[0][2]=2/x*((cx-ab)*da+(x2-b2)*db+(ax-bc)*dc);
        r[1][2]=2/x*((a2-x2)*da+(cx+ab)*db+(bx+ac)*dc);
        r[2][2]=-4*(a*da+b*db);
        return r;
    }
    private static float[][][] shapes(int scale,int k) {
        int nx=FULL_TEMPLATE[0]/scale,ny=FULL_TEMPLATE[1]/scale,nz=FULL_TEMPLATE[2]/scale;
        float[][][] out=new float[nx][ny][nz];
        float divisor=scale*scale*scale;
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) {
            float sum=0.0f;
            for(int dx=0;dx<scale;dx++) for(int dy=0;dy<scale;dy++) for(int dz=0;dz<scale;dz++)
                sum+=(x*scale+dx)>=8?0.0f:((2*(x*scale+dx)+5*(y*scale+dy)+3*(z*scale+dz)+4*k)%11)/10.0f;
            out[x][y][z]=sum/divisor;
        }
        return out;
    }
    private static float[][][][] memberships(int scale) {
        int nx=FULL_IMAGE[0]/scale,ny=FULL_IMAGE[1]/scale,nz=FULL_IMAGE[2]/scale;
        float[][][][] out=new float[nx][ny][nz][3]; float divisor=scale*scale*scale;
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) for(int k=0;k<3;k++) {
            float sum=0.0f;
            for(int dx=0;dx<scale;dx++) for(int dy=0;dy<scale;dy++) for(int dz=0;dz<scale;dz++)
                sum+=((7*(x*scale+dx)+5*(y*scale+dy)+3*(z*scale+dz)+2*k)%17)/16.0f;
            out[x][y][z][k]=sum/divisor;
        }
        return out;
    }
    private static int bounded(int value,int low,int high) { return Math.max(low,Math.min(high,value)); }
    private static float closest(float[][][] a,float x,float y,float z) {
        int nx=a.length,ny=a[0].length,nz=a[0][0].length;
        int x0=bounded((int)Math.floor(x),0,nx-2),y0=bounded((int)Math.floor(y),0,ny-2),z0=bounded((int)Math.floor(z),0,nz-2);
        float alpha=boundedFloat(x-x0,0,1),beta=boundedFloat(y-y0,0,1),gamma=boundedFloat(z-z0,0,1);
        float na=1-alpha,nb=1-beta,ng=1-gamma;
        return na*nb*ng*a[x0][y0][z0]+alpha*nb*ng*a[x0+1][y0][z0]
            +na*beta*ng*a[x0][y0+1][z0]+na*nb*gamma*a[x0][y0][z0+1]
            +alpha*beta*ng*a[x0+1][y0+1][z0]+na*beta*gamma*a[x0][y0+1][z0+1]
            +alpha*nb*gamma*a[x0+1][y0][z0+1]+alpha*beta*gamma*a[x0+1][y0+1][z0+1];
    }
    private static float boundedFloat(float value,float low,float high) { return Math.max(low,Math.min(high,value)); }
    private static float[] priorDerivatives(float[][][] a,float x,float y,float z) {
        int nx=a.length,ny=a[0].length,nz=a[0][0].length;
        if(x<0||x>nx-2||y<0||y>ny-2||z<0||z>nz-2) return new float[3];
        int x0=(int)Math.floor(x),y0=(int)Math.floor(y),z0=(int)Math.floor(z);
        float alpha=x-x0,beta=y-y0,gamma=z-z0,na=1-alpha,nb=1-beta,ng=1-gamma;
        float dx=nb*ng*(a[x0+1][y0][z0]-a[x0][y0][z0])
            +beta*ng*(a[x0+1][y0+1][z0]-a[x0][y0+1][z0])
            +nb*gamma*(a[x0+1][y0][z0+1]-a[x0][y0][z0+1])
            +beta*gamma*(a[x0+1][y0+1][z0+1]-a[x0][y0+1][z0+1]);
        float dy=na*ng*(a[x0][y0+1][z0]-a[x0][y0][z0])
            +alpha*ng*(a[x0+1][y0+1][z0]-a[x0+1][y0][z0])
            +na*gamma*(a[x0][y0+1][z0+1]-a[x0][y0][z0+1])
            +alpha*gamma*(a[x0+1][y0+1][z0+1]-a[x0+1][y0][z0+1]);
        float dz=na*nb*(a[x0][y0][z0+1]-a[x0][y0][z0])
            +alpha*nb*(a[x0+1][y0][z0+1]-a[x0+1][y0][z0])
            +na*beta*(a[x0][y0+1][z0+1]-a[x0][y0+1][z0])
            +alpha*beta*(a[x0+1][y0+1][z0+1]-a[x0+1][y0+1][z0]);
        return new float[]{dx,dy,dz};
    }
    private static float[] mapped(int x,int y,int z,float scale,float[][] rot) {
        float[] center={FULL_IMAGE[0]/2.0f,FULL_IMAGE[1]/2.0f,FULL_IMAGE[2]/2.0f};
        float[] templateCenter={FULL_TEMPLATE[0]/2.0f,FULL_TEMPLATE[1]/2.0f,FULL_TEMPLATE[2]/2.0f};
        float[] q={(x*scale-center[0])*IMAGE_RES[0],(y*scale-center[1])*IMAGE_RES[1],(z*scale-center[2])*IMAGE_RES[2]};
        float[] out=new float[3];
        for(int a=0;a<3;a++) out[a]=((rot[a][0]*q[0]+rot[a][1]*q[1]+rot[a][2]*q[2]+TRANSFORM[a+3])/TEMPLATE_RES[a]+templateCenter[a])/scale;
        return out;
    }
    private static float[][] jacobian(int x,int y,int z,float scale,float[][][] dr) {
        float[] center={FULL_IMAGE[0]/2.0f,FULL_IMAGE[1]/2.0f,FULL_IMAGE[2]/2.0f};
        float[] q={(x*scale-center[0])*IMAGE_RES[0],(y*scale-center[1])*IMAGE_RES[1],(z*scale-center[2])*IMAGE_RES[2]};
        float[][] out=new float[3][6];
        for(int a=0;a<3;a++) {
            for(int p=0;p<3;p++) out[a][p]=(dr[p][a][0]*q[0]+dr[p][a][1]*q[1]+dr[p][a][2]*q[2])/TEMPLATE_RES[a]/scale;
            out[a][a+3]=1.0f/TEMPLATE_RES[a]/scale;
        }
        return out;
    }
    private static float[][] evaluate(int scale) {
        int nx=FULL_IMAGE[0]/scale,ny=FULL_IMAGE[1]/scale,nz=FULL_IMAGE[2]/scale;
        float[][][][] atlas=new float[3][][][]; for(int k=1;k<3;k++) atlas[k]=shapes(scale,k);
        float[][][][] mems=memberships(scale);
        float[][] rot=rotation(TRANSFORM[0],TRANSFORM[1],TRANSFORM[2]);
        float[][][] dr={derivative(TRANSFORM[0],TRANSFORM[1],TRANSFORM[2],1,0,0),derivative(TRANSFORM[0],TRANSFORM[1],TRANSFORM[2],0,1,0),derivative(TRANSFORM[0],TRANSFORM[1],TRANSFORM[2],0,0,1)};
        float cost=0,norm=0; int positive=0,zero=0;
        float[] grad=new float[6],hess=new float[6],ngrad=new float[6],nhess=new float[6];
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) for(int k=1;k<3;k++) {
            float[] p=mapped(x,y,z,scale,rot); float prior=closest(atlas[k],p[0],p[1],p[2]);
            if(prior>0) {
                positive++;
                float[] pd=priorDerivatives(atlas[k],p[0],p[1],p[2]); float[][] j=jacobian(x,y,z,scale,dr); float[] dp=new float[6];
                for(int t=0;t<6;t++) dp[t]=pd[0]*j[0][t]+pd[1]*j[1][t]+pd[2]*j[2][t];
                float m=mems[x][y][z][k];
                cost+=m*m*prior*prior; norm+=prior*prior;
                for(int t=0;t<6;t++) {
                    grad[t]+=2.0f*m*m*prior*dp[t]; hess[t]+=2.0f*m*m*dp[t]*dp[t];
                    ngrad[t]+=2.0f*prior*dp[t]; nhess[t]+=2.0f*dp[t]*dp[t];
                }
            } else zero++;
        }
        return new float[][]{{scale,cost/norm,cost,norm},grad,hess,ngrad,nhess,{positive,zero}};
    }
    private static void append(String key,float[][] values,StringBuilder out) {
        if(out.length()>1) out.append(','); out.append('"').append(key).append("\":{");
        String[] names={"scale_energy_cost_norm","gradient","hessian","norm_gradient","norm_hessian","positive_zero_counts"};
        for(int r=0;r<values.length;r++) { if(r>0)out.append(','); out.append('"').append(names[r]).append("\":[");
            for(int i=0;i<values[r].length;i++){if(i>0)out.append(',');out.append(Float.toString(values[r][i]));} out.append(']'); }
        out.append('}');
    }
    public static void main(String[] args) { StringBuilder out=new StringBuilder("{"); append("scale_2",evaluate(2),out); append("scale_4",evaluate(4),out); System.out.println(out.append('}')); }
}
