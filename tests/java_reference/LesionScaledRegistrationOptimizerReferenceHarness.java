import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Selected-method formula excerpt for LesionToadAtlas's fixed-scale outer loop. */
public final class LesionScaledRegistrationOptimizerReferenceHarness {
    private static final int[] DIMS={8,8,8};
    private static final float[] ZERO={0,0,0,0,0,0};
    private static final int STEP=2, MAX_REJECTIONS=10;
    private static final float FACTOR=1.5f;

    private static float[][][] shape(int clazz) {
        float[][][] a=new float[8][8][8];
        for(int x=0;x<8;x++)for(int y=0;y<8;y++)for(int z=0;z<8;z++) {
            if(clazz==0)a[x][y][z]=0.9f;
            else {int r=(9*x+17*y+19*z+6*x*y+3*y*z)%7;a[x][y][z]=r<3?0.95f:0.05f;}
        }
        return a;
    }
    private static float membership(int x,int y,int z,int clazz) {
        if(clazz==0)return 0.5f;
        int r=(3*x+12*y+13*z+15*x*z+2*x*y)%9;return r<4?0.9f:0.1f;
    }
    private static String vector(float[] a) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}
        return s.append(']').toString();
    }
    private static float[] coefficients(float[] transform,boolean withDerivatives) {
        int n=8; float[] center={4,4,4};
        ParametricTransform model=new ParametricTransform("rigid",4,4,4,1,1,1,8,8,8,4,4,4,1,1,1,8,8,8);
        RotationMatrix rotation=model.computeRotationMatrix(transform);
        float[][] rot=rotation.getMatrix(),da=rotation.derivatives(1,0,0),db=rotation.derivatives(0,1,0),dc=rotation.derivatives(0,0,1);
        float[][][][] maps={shape(0),shape(1)};
        float cost=0,norm=0;float[] grad=new float[6],hess=new float[6];
        for(int x=0;x<n;x+=STEP)for(int y=0;y<n;y+=STEP)for(int z=0;z<n;z+=STEP)for(int k=0;k<2;k++)if(k!=0){
            float[] xi=model.imageToTemplate(x,y,z,transform,rot,1.0f);
            float p=ImageFunctions.linearInterpolation(maps[k],0.0f,xi[0],xi[1],xi[2],8,8,8);
            if(p>0.0f){
                float m=membership(x,y,z,k);
                if(withDerivatives){
                    float dx=ImageFunctions.linearInterpolationXderivative(maps[k],xi[0],xi[1],xi[2],8,8,8);
                    float dy=ImageFunctions.linearInterpolationYderivative(maps[k],xi[0],xi[1],xi[2],8,8,8);
                    float dz=ImageFunctions.linearInterpolationZderivative(maps[k],xi[0],xi[1],xi[2],8,8,8);
                    float[][] dxi=model.imageToTemplateDerivatives(x,y,z,transform,rot,da,db,dc,1.0f);
                    for(int i=0;i<6;i++){float dp=dx*dxi[0][i]+dy*dxi[1][i]+dz*dxi[2][i];hess[i]+=2.0f*m*m*dp*dp;grad[i]+=2.0f*m*m*p*dp;}
                    cost+=m*m*p*p;
                } else cost+=m*m*p*p;
                norm+=p*p;
            }
        }
        float[] out=new float[15];out[0]=cost;out[1]=norm;out[2]=cost/norm;
        for(int i=0;i<6;i++){out[3+i]=grad[i]/norm;out[9+i]=hess[i]/norm;}
        return out;
    }
    private static final class Step {
        float[] transform;float energy,lambda,difference;boolean accepted;int rejected;
        String json;
    }
    private static Step step(float[] current,float lambda) {
        float[] c=coefficients(current,true);float start=c[2],energy=start;int rejected=0;boolean accepted=false;
        StringBuilder trials=new StringBuilder("[");float[] finalTransform=current;
        while(true){
            float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=current[i]+lambda/c[9+i]*c[3+i];
            energy=coefficients(trial,false)[2];if(trials.length()>1)trials.append(',');
            trials.append("{\"transform\":").append(vector(trial)).append(",\"energy\":").append(Float.toString(energy)).append('}');
            if(energy>start){finalTransform=trial;lambda*=FACTOR;accepted=true;break;}
            lambda/=FACTOR;rejected++;if(rejected>MAX_REJECTIONS)break;
        }
        float change=(energy-start)/energy;Step out=new Step();out.transform=finalTransform;out.energy=energy;out.lambda=lambda;out.difference=change;out.accepted=accepted;out.rejected=rejected;
        out.json="{\"trials\":"+trials.append(']')+",\"transform\":"+vector(finalTransform)+",\"energy\":"+Float.toString(energy)+",\"lambda\":"+Float.toString(lambda)+",\"relative_change\":"+Float.toString(change)+",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}';
        return out;
    }
    private static String optimize(String name,int maxIterations,float minChange,float minLambda) {
        float[] trans=ZERO.clone();float lambda=1.0f,diff=1.0f,energy=coefficients(trans,false)[2];int n=0;
        StringBuilder history=new StringBuilder("[");
        while(n<maxIterations && diff>minChange && lambda>minLambda){
            Step s=step(trans,lambda);if(n>0)history.append(',');history.append(s.json);
            trans=s.transform;energy=s.energy;lambda=s.lambda;diff=s.difference;n++;
        }
        String reason=n>=maxIterations?"iteration_limit":(diff<=minChange?"min_energy_change":"lambda_floor");
        return "{\"name\":\""+name+"\",\"max_iterations\":"+maxIterations+",\"min_energy_change\":"+Float.toString(minChange)+",\"min_lambda\":"+Float.toString(minLambda)+",\"iterations\":"+history.append(']')+",\"transform\":"+vector(trans)+",\"energy\":"+Float.toString(energy)+",\"lambda\":"+Float.toString(lambda)+",\"stop_reason\":\""+reason+"\"}";
    }
    public static void main(String[] args) {
        System.out.println("{\"cases\":["+optimize("source_stops",20,1e-6f,0.001f)+","+optimize("iteration_limit",3,1e-6f,0.001f)+","+optimize("minimum_change",20,0.1f,0.001f)+"]}");
    }
}
