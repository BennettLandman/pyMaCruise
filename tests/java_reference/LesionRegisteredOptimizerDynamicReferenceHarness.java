import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/** Selected LesionToadAtlas.registerAllShapes optimizer on supplied state. */
public final class LesionRegisteredOptimizerDynamicReferenceHarness {
    private int nx, ny, nz, ntx, nty, ntz, classes, maxIterations;
    private float[] imageResolution = new float[3], templateResolution = new float[3], transform = new float[6];
    private float[][][][] memberships, shapes;
    private ParametricTransform model;

    private void read(String path) throws Exception {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path)))) {
            nx=in.readInt(); ny=in.readInt(); nz=in.readInt();
            ntx=in.readInt(); nty=in.readInt(); ntz=in.readInt(); classes=in.readInt(); maxIterations=in.readInt();
            for(int i=0;i<3;i++) imageResolution[i]=in.readFloat();
            for(int i=0;i<3;i++) templateResolution[i]=in.readFloat();
            for(int i=0;i<6;i++) transform[i]=in.readFloat();
            memberships=new float[nx][ny][nz][classes];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)
                memberships[x][y][z][k]=in.readFloat();
            shapes=new float[classes][ntx][nty][ntz];
            for(int k=0;k<classes;k++)for(int x=0;x<ntx;x++)for(int y=0;y<nty;y++)for(int z=0;z<ntz;z++)
                shapes[k][x][y][z]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("unexpected trailing input");
        }
        model=new ParametricTransform("rigid",nx/2.0f,ny/2.0f,nz/2.0f,
            imageResolution[0],imageResolution[1],imageResolution[2],nx,ny,nz,
            ntx/2.0f,nty/2.0f,ntz/2.0f,
            templateResolution[0],templateResolution[1],templateResolution[2],ntx,nty,ntz);
    }

    private float[] coefficients(float[] trans, boolean derivative) {
        float[] h=new float[6], g=new float[6]; float cost=0.0f, norm=0.0f;
        float[][] rot=ParametricTransform.computeRotation(trans);
        RotationMatrix rm=new RotationMatrix(); rm.setParameters(trans[0],trans[1],trans[2]);
        float[][] da=rm.derivatives(1.0f,0.0f,0.0f), db=rm.derivatives(0.0f,1.0f,0.0f), dc=rm.derivatives(0.0f,0.0f,1.0f);
        for(int x=0;x<nx;x+=2)for(int y=0;y<ny;y+=2)for(int z=0;z<nz;z+=2)for(int k=1;k<classes;k++) {
            float[] xi=model.imageToTemplate(x,y,z,trans,rot,1.0f);
            float p=ImageFunctions.linearInterpolation(shapes[k],0.0f,xi[0],xi[1],xi[2],ntx,nty,ntz);
            if(p>0.0f) {
                float m=memberships[x][y][z][k];
                if(derivative) {
                    float dx=ImageFunctions.linearInterpolationXderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float dy=ImageFunctions.linearInterpolationYderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float dz=ImageFunctions.linearInterpolationZderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float[][] dxi=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,1.0f);
                    for(int i=0;i<6;i++) {
                        float dp=dx*dxi[0][i]+dy*dxi[1][i]+dz*dxi[2][i];
                        h[i]+=2.0f*m*m*dp*dp; g[i]+=2.0f*m*m*p*dp;
                    }
                    cost+=m*m*p*p;
                } else cost+=m*m*p*p;
                norm+=p*p;
            }
        }
        float[] out=new float[14]; out[0]=cost/norm; out[1]=norm;
        for(int i=0;i<6;i++){out[2+i]=g[i]/norm;out[8+i]=h[i]/norm;}
        return out;
    }

    private float optimize() {
        float lambda=1.0f, difference=1.0f, energy=coefficients(transform,false)[0]; int iteration=0;
        while(iteration<maxIterations && difference>1e-6f && lambda>0.001f) {
            float[] c=coefficients(transform,true); float start=c[0], trialEnergy=start;
            int rejected=0; boolean accepted=false; float[] acceptedTransform=transform;
            while(true) {
                float[] trial=new float[6];
                for(int i=0;i<6;i++)trial[i]=transform[i]+lambda/c[8+i]*c[2+i];
                trialEnergy=coefficients(trial,false)[0];
                if(trialEnergy>start){acceptedTransform=trial;lambda*=1.5f;accepted=true;break;}
                lambda/=1.5f;rejected++;if(rejected>10)break;
            }
            difference=(trialEnergy-start)/trialEnergy;
            if(accepted)transform=acceptedTransform;
            energy=trialEnergy;iteration++;
        }
        return energy;
    }

    private void write(String path, float energy) throws Exception {
        float[][] rotation=ParametricTransform.computeRotation(transform);
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))) {
            out.writeInt(maxIterations); out.writeFloat(energy);
            for(float value:transform)out.writeFloat(value);
            for(float[] row:rotation)for(float value:row)out.writeFloat(value);
            out.writeInt(classes);
            for(int k=0;k<classes;k++) {
                int[] lower={nx,ny,nz}, upper={0,0,0}; int count=0;
                for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
                    float[] xi=model.imageToTemplate(x,y,z,transform,rotation,1.0f);
                    if(ImageFunctions.linearInterpolation(shapes[k],0.0f,xi[0],xi[1],xi[2],ntx,nty,ntz)>0.0f) {
                        count++; if(x<lower[0])lower[0]=x;if(y<lower[1])lower[1]=y;if(z<lower[2])lower[2]=z;
                        if(x>upper[0])upper[0]=x;if(y>upper[1])upper[1]=y;if(z>upper[2])upper[2]=z;
                    }
                }
                for(int value:lower)out.writeInt(value);
                for(int value:upper)out.writeInt(value);
                out.writeInt(count);
            }
        }
    }

    public static void main(String[] args) throws Exception {
        if(args.length!=2)throw new IllegalArgumentException("input.bin output.bin required");
        LesionRegisteredOptimizerDynamicReferenceHarness h=new LesionRegisteredOptimizerDynamicReferenceHarness();
        h.read(args[0]); h.write(args[1],h.optimize());
    }
}
