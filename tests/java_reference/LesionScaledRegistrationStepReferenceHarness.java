import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Selected-method formula excerpt for one scaled LesionToadAtlas update. */
public final class LesionScaledRegistrationStepReferenceHarness {
    private static final int[] IMAGE_DIMS={12,12,12}, TEMPLATE_DIMS={12,12,12};
    private static final float[] IMAGE_RES={1.0f,1.0f,1.0f}, TEMPLATE_RES={1.0f,1.0f,1.0f};
    private static final float[] INITIAL={0,0,0,0,0,0};
    private static final int CLASSES=2, SUBSAMPLE=2, MAX_REJECTIONS=10;
    private static final float LAMBDA_FACTOR=1.5f;

    private static int[] coarse(int[] dims,int scale) { return new int[]{dims[0]/scale,dims[1]/scale,dims[2]/scale}; }
    private static float[][][] shape(int k,int[] dims) {
        float[][][] out=new float[dims[0]][dims[1]][dims[2]];
        for(int x=0;x<dims[0];x++) for(int y=0;y<dims[1];y++) for(int z=0;z<dims[2];z++) {
            if(k==0) out[x][y][z]=0.9f;
            else {
                int r=(15*x+4*y+10*z+5*x*y+8*y*z+3*x*z+3)%17;
                out[x][y][z]=0.05f+r/20.0f;
            }
        }
        return out;
    }
    private static float membership(int x,int y,int z,int k,int scale) {
        if(k==0) return 0.5f;
        int a=scale==1?23:13,b=scale==1?15:7,c=scale==1?8:5,d=scale==1?14:3;
        int r=(a*x+b*y+c*z+d)%23;
        return 0.04f+r/25.0f;
    }
    private static String vector(float[] a) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}
        return s.append(']').toString();
    }
    private static float[] compute(float[] trans,int scale,boolean coefficients) {
        int[] in=coarse(IMAGE_DIMS,scale), tpl=coarse(TEMPLATE_DIMS,scale);
        float[] ic={IMAGE_DIMS[0]/2.0f,IMAGE_DIMS[1]/2.0f,IMAGE_DIMS[2]/2.0f};
        float[] tc={TEMPLATE_DIMS[0]/2.0f,TEMPLATE_DIMS[1]/2.0f,TEMPLATE_DIMS[2]/2.0f};
        ParametricTransform model=new ParametricTransform("rigid",ic[0],ic[1],ic[2],
            IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],IMAGE_DIMS[0],IMAGE_DIMS[1],IMAGE_DIMS[2],
            tc[0],tc[1],tc[2],TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],
            TEMPLATE_DIMS[0],TEMPLATE_DIMS[1],TEMPLATE_DIMS[2]);
        RotationMatrix rotation=model.computeRotationMatrix(trans);
        float[][] rot=rotation.getMatrix();
        float[][] da=rotation.derivatives(1,0,0), db=rotation.derivatives(0,1,0), dc=rotation.derivatives(0,0,1);
        float[][][][] maps=new float[CLASSES][][][];
        for(int k=0;k<CLASSES;k++) maps[k]=shape(k,tpl);
        float cost=0,norm=0; float[] grad=new float[6],hess=new float[6];
        for(int x=0;x<in[0];x+=SUBSAMPLE) for(int y=0;y<in[1];y+=SUBSAMPLE) for(int z=0;z<in[2];z+=SUBSAMPLE) {
            for(int k=0;k<CLASSES;k++) if(k!=0) {
                float[] xi=model.imageToTemplate(x,y,z,trans,rot,(float)scale);
                float p=ImageFunctions.linearInterpolation(maps[k],0.0f,xi[0],xi[1],xi[2],tpl[0],tpl[1],tpl[2]);
                if(p>0.0f) {
                    if(coefficients) {
                        float dx=ImageFunctions.linearInterpolationXderivative(maps[k],xi[0],xi[1],xi[2],tpl[0],tpl[1],tpl[2]);
                        float dy=ImageFunctions.linearInterpolationYderivative(maps[k],xi[0],xi[1],xi[2],tpl[0],tpl[1],tpl[2]);
                        float dz=ImageFunctions.linearInterpolationZderivative(maps[k],xi[0],xi[1],xi[2],tpl[0],tpl[1],tpl[2]);
                        float[][] dxi=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,(float)scale);
                        float m=membership(x,y,z,k,scale);
                        for(int i=0;i<6;i++) {
                            float dp=dx*dxi[0][i]+dy*dxi[1][i]+dz*dxi[2][i];
                            hess[i]+=2.0f*m*m*dp*dp;
                            grad[i]+=2.0f*m*m*p*dp;
                        }
                        cost+=m*m*p*p;
                    }
                    else cost+=membership(x,y,z,k,scale)*membership(x,y,z,k,scale)*p*p;
                    norm+=p*p;
                }
            }
        }
        float[] out=new float[15]; out[0]=cost;out[1]=norm;out[2]=cost/norm;
        for(int i=0;i<6;i++){out[3+i]=grad[i]/norm;out[9+i]=hess[i]/norm;}
        return out;
    }
    private static String run(int scale,float startLambda) {
        float[] current=INITIAL.clone();
        float[] initial=compute(current,scale,true); float e0=initial[2],energy=e0,lambda=startLambda;
        StringBuilder trials=new StringBuilder("["); int rejected=0; boolean accepted=false;
        while(true) {
            float[] trial=new float[6];
            for(int i=0;i<6;i++) trial[i]=current[i]+lambda/initial[9+i]*initial[3+i];
            energy=compute(trial,scale,false)[2];
            if(trials.length()>1)trials.append(',');
            trials.append("{\"transform\":").append(vector(trial)).append(",\"energy\":").append(Float.toString(energy)).append('}');
            if(energy>e0){current=trial;lambda=lambda*LAMBDA_FACTOR;accepted=true;break;}
            lambda=lambda/LAMBDA_FACTOR;rejected++;
            if(rejected>MAX_REJECTIONS)break;
        }
        float relative=(energy-e0)/energy;
        return "{\"scale\":"+scale+",\"initial_lambda\":"+Float.toString(startLambda)
            +",\"initial_energy\":"+Float.toString(e0)+",\"trials\":"+trials.append(']')
            +",\"transform\":"+vector(current)+",\"energy\":"+Float.toString(energy)
            +",\"lambda\":"+Float.toString(lambda)+",\"relative_change\":"+Float.toString(relative)
            +",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}';
    }
    public static void main(String[] args) {
        System.out.println("{\"cases\":["+run(1,0.000001f)+","+run(2,0.0001f)+"]}");
    }
}
