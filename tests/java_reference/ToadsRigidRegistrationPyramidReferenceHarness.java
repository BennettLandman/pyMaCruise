import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Selected-source TOADS coarse-to-fine rigid registration reference. */
public final class ToadsRigidRegistrationPyramidReferenceHarness {
    private static final int MAX_REJECTIONS=10;
    private static final float FACTOR=1.5f;
    private static final float[] IMAGE_RES={0.8f,1.2f,1.5f}, TEMPLATE_RES={0.9f,1.1f,1.3f};
    private static final float[] INITIAL={0.04f,-0.025f,0.015f,0.03f,-0.02f,0.01f};

    private static final class Volume {
        final int nx,ny,nz,nc; final float[][][][] data;
        Volume(Path file) throws IOException {
            byte[] bytes;
            try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(file)); ByteArrayOutputStream out=new ByteArrayOutputStream()) {
                in.transferTo(out); bytes=out.toByteArray();
            }
            ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
            int dims=h.getShort(40), datatype=h.getShort(70); float offset=h.getFloat(108);
            nx=h.getShort(42); ny=h.getShort(44); nz=h.getShort(46); nc=dims>3?h.getShort(48):1;
            if(datatype!=16 || nx<1 || ny<1 || nz<1 || nc<1) throw new IOException("expected float32 NIfTI-1");
            data=new float[nx][ny][nz][nc]; int start=(int)offset;
            for(int c=0;c<nc;c++) for(int z=0;z<nz;z++) for(int y=0;y<ny;y++) for(int x=0;x<nx;x++) {
                int i=x+nx*(y+ny*(z+nz*c)); data[x][y][z][c]=h.getFloat(start+4*i);
            }
        }
    }
    private static String vector(float[] a) {
        StringBuilder s=new StringBuilder("["); for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));} return s.append(']').toString();
    }
    private static float[] measure(float[] trans,float[][][][] mem,float[][][][] shapes,int scale,boolean derivatives,int[] fullImage,int[] fullTemplate) {
        int nx=mem.length,ny=mem[0].length,nz=mem[0][0].length;
        ParametricTransform model=new ParametricTransform("rigid",fullImage[0]/2.0f,fullImage[1]/2.0f,fullImage[2]/2.0f,
            IMAGE_RES[0],IMAGE_RES[1],IMAGE_RES[2],fullImage[0],fullImage[1],fullImage[2],
            fullTemplate[0]/2.0f,fullTemplate[1]/2.0f,fullTemplate[2]/2.0f,
            TEMPLATE_RES[0],TEMPLATE_RES[1],TEMPLATE_RES[2],fullTemplate[0],fullTemplate[1],fullTemplate[2]);
        RotationMatrix r=model.computeRotationMatrix(trans); float[][] rot=r.getMatrix();
        // MultipleToadAtlas initializes precompute=true for a linear transform
        // and computeExactScaledRegistrationCoefficients refreshes this matrix
        // at the active pyramid scale before sampling each voxel.
        float[][] matrix=new float[3][4];
        model.precomputeImageToTemplateMatrix(matrix,trans,rot,(float)scale);
        float[][] da=r.derivatives(1,0,0),db=r.derivatives(0,1,0),dc=r.derivatives(0,0,1);
        float cost=0,norm=0; float[] grad=new float[6],hess=new float[6];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=1;k<mem[0][0][0].length;k++) {
            float[] pxyz=new float[3];
            for(int axis=0;axis<3;axis++)
                pxyz[axis]=matrix[axis][0]*x+matrix[axis][1]*y+matrix[axis][2]*z+matrix[axis][3];
            float p=ImageFunctions.linearClosestInterpolation(shapes[k],pxyz[0],pxyz[1],pxyz[2],shapes[k].length,shapes[k][0].length,shapes[k][0][0].length);
            if(p<=0)continue;
            float m=mem[x][y][z][k]; float m2=m*m; float p2=p*p;
            // Match MultipleToadDeformableAtlas.registrationCost: Java
            // evaluates m*m*p*p left-to-right, rounding at each float step.
            cost+=m*m*p*p; norm+=p2;
            if(derivatives) {
                float px=ImageFunctions.linearInterpolationXderivative(shapes[k],pxyz[0],pxyz[1],pxyz[2],shapes[k].length,shapes[k][0].length,shapes[k][0][0].length);
                float py=ImageFunctions.linearInterpolationYderivative(shapes[k],pxyz[0],pxyz[1],pxyz[2],shapes[k].length,shapes[k][0].length,shapes[k][0][0].length);
                float pz=ImageFunctions.linearInterpolationZderivative(shapes[k],pxyz[0],pxyz[1],pxyz[2],shapes[k].length,shapes[k][0].length,shapes[k][0][0].length);
                float[][] jac=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,(float)scale);
                for(int i=0;i<6;i++){float dp=px*jac[0][i]+py*jac[1][i]+pz*jac[2][i];grad[i]+=2.0f*m2*p*dp;hess[i]+=2.0f*m2*dp*dp;}
            }
        }
        float[] out=new float[15];out[0]=cost;out[1]=norm;out[2]=cost/norm;
        for(int i=0;i<6;i++){out[3+i]=grad[i];out[9+i]=hess[i];}
        return out;
    }
    private static final class Step {
        float[] transform; float energy,lambda,difference; boolean accepted; int rejected;
        String json;
    }
    private static Step update(float[] current,float lambda,float[][][][] mem,float[][][][] shapes,int scale,int[] fi,int[] ft) {
        float[] c=measure(current,mem,shapes,scale,true,fi,ft); float start=c[2],energy=start;
        float[] gradient=new float[6],hessian=new float[6];
        for(int i=0;i<6;i++){hessian[i]=c[9+i]/c[1];gradient[i]=hessian[i]/c[1];}
        int rejected=0; boolean accepted=false; float[] finalTransform=current; StringBuilder trials=new StringBuilder("[");
        while(true) {
            float[] trial=new float[6]; for(int i=0;i<6;i++)trial[i]=current[i]+lambda/hessian[i]*gradient[i];
            energy=measure(trial,mem,shapes,scale,false,fi,ft)[2]; if(trials.length()>1)trials.append(',');
            trials.append("{\"transform\":").append(vector(trial)).append(",\"energy\":").append(Float.toString(energy)).append('}');
            if(energy>start){finalTransform=trial;lambda*=FACTOR;accepted=true;break;}
            lambda/=FACTOR;rejected++;if(rejected>MAX_REJECTIONS)break;
        }
        float difference=(energy-start)/energy; Step out=new Step();out.transform=finalTransform;out.energy=energy;out.lambda=lambda;out.difference=difference;out.accepted=accepted;out.rejected=rejected;
        out.json="{\"trials\":"+trials.append(']')+",\"transform\":"+vector(finalTransform)+",\"energy\":"+energy+",\"lambda\":"+lambda+",\"relative_change\":"+difference+",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}'; return out;
    }
    private static String pyramid(Volume image,Volume atlas) {
        float[] trans=INITIAL.clone(); int[] fi={image.nx,image.ny,image.nz},ft={atlas.nx,atlas.ny,atlas.nz};
        StringBuilder stages=new StringBuilder("["); boolean first=true;
        for(int scale:new int[]{4,2}) {
            float[][][][] mem=ImageFunctions.subsample(image.data,image.nx,image.ny,image.nz,image.nc,scale);
            float[][][][] shp=new float[atlas.nc][][][]; for(int k=0;k<atlas.nc;k++){float[][][] channel=new float[atlas.nx][atlas.ny][atlas.nz];for(int x=0;x<atlas.nx;x++)for(int y=0;y<atlas.ny;y++)for(int z=0;z<atlas.nz;z++)channel[x][y][z]=atlas.data[x][y][z][k];shp[k]=ImageFunctions.subsample(channel,atlas.nx,atlas.ny,atlas.nz,scale);}
            float[] levelStart=trans.clone(); float lambda=1.0f,difference=1.0f,energy=measure(trans,mem,shp,scale,false,fi,ft)[2];
            List<String> iterations=new ArrayList<>();
            for(int n=0;n<2 && difference>1e-6f && lambda>0.001f;n++){Step s=update(trans,lambda,mem,shp,scale,fi,ft);iterations.add(s.json);trans=s.transform;energy=s.energy;lambda=s.lambda;difference=s.difference;}
            if(!first)stages.append(',');first=false;
            stages.append("{\"scale\":").append(scale).append(",\"input_dims\":[").append(mem.length).append(',').append(mem[0].length).append(',').append(mem[0][0].length)
                .append("],\"template_dims\":[").append(shp[1].length).append(',').append(shp[1][0].length).append(',').append(shp[1][0][0].length)
                .append("],\"initial_transform\":").append(vector(levelStart)).append(",\"iterations\":[").append(String.join(",",iterations))
                .append("],\"transform\":").append(vector(trans)).append(",\"energy\":").append(energy).append(",\"lambda\":").append(lambda).append(",\"relative_change\":").append(difference).append('}');
        }
        return "{\"initial_transform\":"+vector(INITIAL)+",\"stages\":"+stages.append(']')+",\"final_transform\":"+vector(trans)+'}';
    }
    public static void main(String[] args) throws Exception {
        Volume mem=new Volume(Path.of(args[0])),atlas=new Volume(Path.of(args[1]));
        System.out.println(pyramid(mem,atlas));
    }
}
