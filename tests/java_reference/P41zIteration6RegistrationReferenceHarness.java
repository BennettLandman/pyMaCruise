import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;

/**
 * Selected LesionToadAtlas.registerAllShapes() formula excerpt for the fixed
 * P41z iteration-6 state. Input and output arrays are in RAS XYZ order.
 * The harness uses supplied JIST ParametricTransform and ImageFunctions.
 */
public final class P41zIteration6RegistrationReferenceHarness {
    private static final int SUBSAMPLE=2, IT_MAX=10, IT_SUPP=10;
    private static final float MIN_EDIFF=1.0e-6f, MIN_LAMBDA=0.001f, FACTOR=1.5f;
    private int nx,ny,nz,tx,ty,tz,classes;
    private float[] imageRes,templateRes,initial;
    private float[][][][] shapes,memberships;
    private ParametricTransform model;

    private void read(String path) throws Exception {
        try(DataInputStream in=new DataInputStream(new BufferedInputStream(new FileInputStream(path)))) {
            nx=in.readInt();ny=in.readInt();nz=in.readInt();tx=in.readInt();ty=in.readInt();tz=in.readInt();classes=in.readInt();
            if(nx<3||ny<3||nz<3||tx<2||ty<2||tz<2||classes<2)throw new IllegalArgumentException("invalid dimensions");
            imageRes=new float[3];templateRes=new float[3];initial=new float[6];
            for(int i=0;i<3;i++)imageRes[i]=in.readFloat();
            for(int i=0;i<3;i++)templateRes[i]=in.readFloat();
            for(int i=0;i<6;i++)initial[i]=in.readFloat();
            shapes=new float[classes][tx][ty][tz];
            memberships=new float[nx][ny][nz][classes];
            for(int k=0;k<classes;k++)for(int x=0;x<tx;x++)for(int y=0;y<ty;y++)for(int z=0;z<tz;z++)shapes[k][x][y][z]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)memberships[x][y][z][k]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("trailing input bytes");
        }
        model=new ParametricTransform("rigid",nx/2.0f,ny/2.0f,nz/2.0f,imageRes[0],imageRes[1],imageRes[2],nx,ny,nz,
            tx/2.0f,ty/2.0f,tz/2.0f,templateRes[0],templateRes[1],templateRes[2],tx,ty,tz);
    }

    private float[] coefficients(float[] transform,boolean derivative) {
        float[] gradient=new float[6],hessian=new float[6];float cost=0.0f,norm=0.0f;
        float[][] rotation=ParametricTransform.computeRotation(transform);
        RotationMatrix rm=new RotationMatrix();rm.setParameters(transform[0],transform[1],transform[2]);
        float[][] da=rm.derivatives(1,0,0),db=rm.derivatives(0,1,0),dc=rm.derivatives(0,0,1);
        for(int x=0;x<nx;x+=SUBSAMPLE)for(int y=0;y<ny;y+=SUBSAMPLE)for(int z=0;z<nz;z+=SUBSAMPLE)for(int k=1;k<classes;k++) {
            float[] point=model.imageToTemplate(x,y,z,transform,rotation,1.0f);
            float prior=ImageFunctions.linearInterpolation(shapes[k],0.0f,point[0],point[1],point[2],tx,ty,tz);
            if(prior>0.0f) {
                float membership=memberships[x][y][z][k];
                if(derivative) {
                    float dx=ImageFunctions.linearInterpolationXderivative(shapes[k],point[0],point[1],point[2],tx,ty,tz);
                    float dy=ImageFunctions.linearInterpolationYderivative(shapes[k],point[0],point[1],point[2],tx,ty,tz);
                    float dz=ImageFunctions.linearInterpolationZderivative(shapes[k],point[0],point[1],point[2],tx,ty,tz);
                    float[][] dpoint=model.imageToTemplateDerivatives(x,y,z,transform,rotation,da,db,dc,1.0f);
                    for(int i=0;i<6;i++) {
                        float dp=dx*dpoint[0][i]+dy*dpoint[1][i]+dz*dpoint[2][i];
                        hessian[i]+=2.0f*membership*membership*dp*dp;
                        gradient[i]+=2.0f*membership*membership*prior*dp;
                    }
                    cost+=membership*membership*prior*prior;
                } else cost+=membership*membership*prior*prior;
                norm+=prior*prior;
            }
        }
        float[] out=new float[14];out[0]=cost/norm;out[1]=norm;
        for(int i=0;i<6;i++){out[2+i]=gradient[i]/norm;out[8+i]=hessian[i]/norm;}
        return out;
    }

    private static String vector(float[] values) {
        StringBuilder out=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.toString(values[i]));}
        return out.append(']').toString();
    }
    private static String matrix(float[][] values) {
        StringBuilder out=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(vector(values[i]));}
        return out.append(']').toString();
    }
    private String bounds(float[] transform) {
        float[][] rotation=ParametricTransform.computeRotation(transform);StringBuilder out=new StringBuilder("[");
        for(int k=0;k<classes;k++) {
            if(k>0)out.append(',');int[] low={nx,ny,nz},high={0,0,0};int count=0;
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
                float[] point=model.imageToTemplate(x,y,z,transform,rotation,1.0f);
                if(ImageFunctions.linearInterpolation(shapes[k],0.0f,point[0],point[1],point[2],tx,ty,tz)>0.0f) {
                    count++;if(x<low[0])low[0]=x;if(y<low[1])low[1]=y;if(z<low[2])low[2]=z;
                    if(x>high[0])high[0]=x;if(y>high[1])high[1]=y;if(z>high[2])high[2]=z;
                }
            }
            out.append("{\"lower_xyz\":").append(vector(new float[]{low[0],low[1],low[2]}))
               .append(",\"upper_xyz\":").append(vector(new float[]{high[0],high[1],high[2]}))
               .append(",\"support_count\":").append(count).append('}');
        }
        return out.append(']').toString();
    }
    private String run() {
        float[] transform=initial.clone();float lambda=1.0f,diff=1.0f,energy=coefficients(transform,false)[0];
        StringBuilder steps=new StringBuilder("[");int iteration=0;
        while(iteration<IT_MAX&&diff>MIN_EDIFF&&lambda>MIN_LAMBDA) {
            if(iteration>0)steps.append(',');float[] coefficient=coefficients(transform,true);float start=coefficient[0],trialEnergy=start;
            int rejected=0;boolean accepted=false;float[] acceptedTransform=transform;
            StringBuilder trials=new StringBuilder("[");
            while(true) {
                float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=transform[i]+lambda/coefficient[8+i]*coefficient[2+i];
                trialEnergy=coefficients(trial,false)[0];
                if(trials.length()>1)trials.append(',');trials.append("{\"transform\":").append(vector(trial)).append(",\"energy\":").append(Float.toString(trialEnergy)).append('}');
                if(trialEnergy>start){acceptedTransform=trial;lambda*=FACTOR;accepted=true;break;}
                lambda/=FACTOR;rejected++;if(rejected>IT_SUPP)break;
            }
            diff=(trialEnergy-start)/trialEnergy;if(accepted)transform=acceptedTransform;energy=trialEnergy;
            steps.append("{\"trials\":").append(trials.append(']')).append(",\"accepted\":").append(accepted)
                .append(",\"rejected_trials\":").append(rejected).append(",\"transform\":").append(vector(transform))
                .append(",\"energy\":").append(Float.toString(energy)).append(",\"relative_change\":").append(Float.toString(diff))
                .append(",\"lambda\":").append(Float.toString(lambda)).append('}');iteration++;
        }
        String reason=iteration>=IT_MAX?"iteration_limit":(diff<=MIN_EDIFF?"min_energy_change":"lambda_floor");
        return "{\"iterations\":"+steps.append(']')+",\"transform\":"+vector(transform)+",\"energy\":"+energy
            +",\"lambda\":"+lambda+",\"stop_reason\":\""+reason+"\",\"rotation\":"+matrix(ParametricTransform.computeRotation(transform))
            +",\"bounds\":"+bounds(transform)+"}";
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=1)throw new IllegalArgumentException("input.bin");
        P41zIteration6RegistrationReferenceHarness h=new P41zIteration6RegistrationReferenceHarness();h.read(args[0]);
        System.out.println(h.run());
    }
}
