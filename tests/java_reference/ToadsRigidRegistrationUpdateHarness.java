import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import edu.jhmi.rad.medic.libraries.ImageFunctions;

/** Selected-source exact rigid coefficients and registerGradientDescent excerpt. */
public final class ToadsRigidRegistrationUpdateHarness {
    private static final int NX=4,NY=4,NZ=4,TX=5,TY=5,TZ=5,NC=3,ND=6;
    private static final float[][][][] mem=new float[NX][NY][NZ][NC];
    private static final float[][][][] shape=new float[NC][TX][TY][TZ];
    private static final ParametricTransform model=new ParametricTransform(
        "rigid",NX/2.0f,NY/2.0f,NZ/2.0f,1,1,1,NX,NY,NZ,
        TX/2.0f,TY/2.0f,TZ/2.0f,1,1,1,TX,TY,TZ);

    private static void initialize(){
        for(int k=0;k<NC;k++)for(int x=0;x<TX;x++)for(int y=0;y<TY;y++)for(int z=0;z<TZ;z++)
            shape[k][x][y][z]=0.18f+0.035f*x+0.025f*y+0.02f*z+0.03f*k+0.004f*x*y;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            mem[x][y][z][0]=0.1f;
            mem[x][y][z][1]=0.25f+0.04f*x+0.02f*y;
            mem[x][y][z][2]=0.2f+0.015f*z+0.01f*x;
        }
    }
    private static int bits(float x){return Float.floatToIntBits(x);}
    private static String floatBits(float[] v){StringBuilder b=new StringBuilder("[");for(int i=0;i<v.length;i++){if(i>0)b.append(',');b.append(bits(v[i]));}return b.append(']').toString();}
    private static float[] map(float[][] matrix,int x,int y,int z){
        return new float[]{matrix[0][0]*x+matrix[0][1]*y+matrix[0][2]*z+matrix[0][3],
            matrix[1][0]*x+matrix[1][1]*y+matrix[1][2]*z+matrix[1][3],
            matrix[2][0]*x+matrix[2][1]*y+matrix[2][2]*z+matrix[2][3]};
    }
    private static float energy(float[] trans){
        float[][] rot=model.computeRotation(trans),matrix=new float[3][4];
        model.precomputeImageToTemplateMatrix(matrix,trans,rot,1.0f);
        float cost=0,norm=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=1;k<NC;k++){
            float[] xi=map(matrix,x,y,z);
            float p=ImageFunctions.linearClosestInterpolation(shape[k],xi[0],xi[1],xi[2],TX,TY,TZ);
            if(p>0){cost+=mem[x][y][z][k]*mem[x][y][z][k]*p*p;norm+=p*p;}
        }
        return cost/norm;
    }
    private static float[][][] coefficients(float[] trans){
        float cost=0,norm=0;float[] gradient=new float[ND],hessian=new float[ND],ngrad=new float[ND],nhess=new float[ND];
        RotationMatrix r=model.computeRotationMatrix(trans);float[][] rot=r.getMatrix();
        float[][] dRa=r.derivatives(1,0,0),dRb=r.derivatives(0,1,0),dRc=r.derivatives(0,0,1);
        float[][] matrix=new float[3][4];model.precomputeImageToTemplateMatrix(matrix,trans,rot,1.0f);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=1;k<NC;k++){
            float[] xi=map(matrix,x,y,z);
            float p=ImageFunctions.linearClosestInterpolation(shape[k],xi[0],xi[1],xi[2],TX,TY,TZ);
            if(p<=0)continue;
            float px=ImageFunctions.linearInterpolationXderivative(shape[k],xi[0],xi[1],xi[2],TX,TY,TZ);
            float py=ImageFunctions.linearInterpolationYderivative(shape[k],xi[0],xi[1],xi[2],TX,TY,TZ);
            float pz=ImageFunctions.linearInterpolationZderivative(shape[k],xi[0],xi[1],xi[2],TX,TY,TZ);
            float[][] dx=model.imageToTemplateDerivatives(x,y,z,trans,rot,dRa,dRb,dRc,1.0f);
            float[] dp=new float[ND];
            for(int i=0;i<ND;i++)dp[i]=px*dx[0][i]+py*dx[1][i]+pz*dx[2][i];
            float m=mem[x][y][z][k];cost+=m*m*p*p;norm+=p*p;
            for(int i=0;i<ND;i++){
                gradient[i]+=2.0f*m*m*p*dp[i];hessian[i]+=2.0f*m*m*dp[i]*dp[i];
                ngrad[i]+=2.0f*p*dp[i];nhess[i]+=2.0f*dp[i]*dp[i];
            }
        }
        float[] sourceGradient=new float[ND],sourceHessian=new float[ND];
        for(int i=0;i<ND;i++){
            sourceHessian[i]=hessian[i]/norm;
            // Retain selected source call order: pass already-normalized hessian as gradient arg.
            sourceGradient[i]=sourceHessian[i]/norm;
        }
        return new float[][][]{{{cost,norm}},{gradient,hessian,ngrad,nhess},{sourceGradient,sourceHessian}};
    }
    private static void run(String name,float[] initial,StringBuilder out){
        float[] trans=initial.clone();float[][][] coeff=coefficients(trans);float e0=energy(trans),lambda=1.0f,e=e0;
        List<String> trials=new ArrayList<>();int rejected=0;boolean accepted=false;
        while(true){
            float[] trial=new float[ND];for(int i=0;i<ND;i++)trial[i]=trans[i]+lambda/coeff[2][1][i]*coeff[2][0][i];
            e=energy(trial);boolean improves=e>e0;
            trials.add("["+bits(lambda)+","+floatBits(trial)+","+bits(e)+","+(improves?1:0)+"]");
            if(improves){trans=trial;lambda=lambda*1.5f;accepted=true;break;}
            lambda=lambda/1.5f;rejected++;if(rejected>10)break;
        }
        float relative=(e-e0)/e;
        if(out.charAt(out.length()-1)!='{')out.append(',');
        out.append("\"").append(name).append("\":{\"initial_transform\":").append(floatBits(initial))
            .append(",\"initial_energy\":").append(bits(e0)).append(",\"cost_norm\":").append(floatBits(coeff[0][0]))
            .append(",\"raw_gradient\":").append(floatBits(coeff[1][0])).append(",\"raw_hessian\":").append(floatBits(coeff[1][1]))
            .append(",\"norm_gradient\":").append(floatBits(coeff[1][2])).append(",\"norm_hessian\":").append(floatBits(coeff[1][3]))
            .append(",\"gradient\":").append(floatBits(coeff[2][0])).append(",\"hessian\":").append(floatBits(coeff[2][1]))
            .append(",\"trials\":[").append(String.join(",",trials)).append("]")
            .append(",\"transform\":").append(floatBits(trans)).append(",\"energy\":").append(bits(e))
            .append(",\"lambda\":").append(bits(lambda)).append(",\"relative_change\":").append(bits(relative))
            .append(",\"accepted\":").append(accepted).append(",\"rejected\":").append(rejected).append('}');
    }
    /** Selected registerAllShapes outer loop around registerGradientDescent. */
    private static void runMultiple(StringBuilder out){
        float[] trans={0,0,0,0.05f,-0.03f,0.02f};
        float lambda=1.0f,diff=1.0f;int count=0;
        List<String> updates=new ArrayList<>();
        while(count<8 && diff>0.0f && lambda>0.001f){
            float[][][] coeff=coefficients(trans);float e0=energy(trans),e=e0;
            List<String> trials=new ArrayList<>();int rejected=0;boolean accepted=false;
            while(true){
                float[] trial=new float[ND];
                for(int i=0;i<ND;i++)trial[i]=trans[i]+lambda/coeff[2][1][i]*coeff[2][0][i];
                e=energy(trial);boolean improves=e>e0;
                trials.add("["+bits(lambda)+","+floatBits(trial)+","+bits(e)+","+(improves?1:0)+"]");
                if(improves){trans=trial;lambda=lambda*1.5f;accepted=true;break;}
                lambda=lambda/1.5f;rejected++;if(rejected>10)break;
            }
            diff=(e-e0)/e;count++;
            updates.add("{\"trials\":["+String.join(",",trials)+"],\"transform\":"+floatBits(trans)
                +",\"energy\":"+bits(e)+",\"lambda\":"+bits(lambda)+",\"relative_change\":"+bits(diff)
                +",\"accepted\":"+accepted+",\"rejected\":"+rejected+"}");
        }
        if(out.charAt(out.length()-1)!='{')out.append(',');
        out.append("\"multi_step\":{\"initial_transform\":").append(floatBits(new float[]{0,0,0,0.05f,-0.03f,0.02f}))
            .append(",\"updates\":[").append(String.join(",",updates)).append("]}");
    }
    public static void main(String[] args){
        initialize();StringBuilder out=new StringBuilder("{\"source\":\"VandyCRUISE MultipleToadDeformableAtlas.computeExactRegistrationCoefficients + registerGradientDescent\",\"source_sha256\":\"e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2\",\"parametric_transform_sha256\":\"22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac\",\"rotation_matrix_sha256\":\"144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0\",\"image_functions_sha256\":\"be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892\",\"cases\":{");
        run("accepted",new float[]{0,0,0,0.05f,-0.03f,0.02f},out);
        run("backtracked",new float[]{0.169f,0.064f,0.126f,-0.247f,0.078f,0.288f},out);
        runMultiple(out);
        System.out.println(out.append("}}"));
    }
}
