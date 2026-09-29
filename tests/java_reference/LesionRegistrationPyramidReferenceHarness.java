import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;

/** Selected-method formula excerpt for LesionToadAtlas pyramid handoff. */
public final class LesionRegistrationPyramidReferenceHarness {
    private static final int N=33, CLASSES=2, STEP=2, MAX_REJECTIONS=10;
    private static final float FACTOR=1.5f;
    private static float[][][][] fullShapes() {
        float[][][][] a=new float[CLASSES][N][N][N];
        for(int k=0;k<CLASSES;k++)for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            if(k==0)a[k][x][y][z]=0.9f;
            else {int r=(3*x+5*y+7*z+11*x*y+13*y*z)%7;a[k][x][y][z]=r<3?0.95f:0.05f;}
        }
        return a;
    }
    private static float[][][][] fullMemberships() {
        float[][][][] a=new float[N][N][N][CLASSES];
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            a[x][y][z][0]=0.5f;
            int r=(7*x+11*y+5*z+13*x*z+17*x*y)%9;
            a[x][y][z][1]=r<4?0.9f:0.1f;
        }
        return a;
    }
    private static String vector(float[] a) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}
        return s.append(']').toString();
    }
    private static String jsonNumber(float value) {
        return Float.isFinite(value) ? Float.toString(value) : "null";
    }
    private static float[] compute(float[] trans,float[][][][] mem,float[][][][] shapes,int scale,boolean derivatives) {
        int nx=mem.length,ny=mem[0].length,nz=mem[0][0].length;
        int ntx=shapes[1].length,nty=shapes[1][0].length,ntz=shapes[1][0][0].length;
        float[] center={N/2.0f,N/2.0f,N/2.0f};
        ParametricTransform model=new ParametricTransform("rigid",center[0],center[1],center[2],1,1,1,N,N,N,
            center[0],center[1],center[2],1,1,1,N,N,N);
        RotationMatrix r=model.computeRotationMatrix(trans);float[][] rot=r.getMatrix();
        float[][] da=r.derivatives(1,0,0),db=r.derivatives(0,1,0),dc=r.derivatives(0,0,1);
        float cost=0,norm=0;float[] grad=new float[6],hess=new float[6];
        for(int x=0;x<nx;x+=STEP)for(int y=0;y<ny;y+=STEP)for(int z=0;z<nz;z+=STEP)for(int k=0;k<CLASSES;k++)if(k!=0) {
            float[] xi=model.imageToTemplate(x,y,z,trans,rot,(float)scale);
            float p=ImageFunctions.linearInterpolation(shapes[k],0.0f,xi[0],xi[1],xi[2],ntx,nty,ntz);
            if(p>0.0f) {
                float m=mem[x][y][z][k];
                if(derivatives) {
                    float dx=ImageFunctions.linearInterpolationXderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float dy=ImageFunctions.linearInterpolationYderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float dz=ImageFunctions.linearInterpolationZderivative(shapes[k],xi[0],xi[1],xi[2],ntx,nty,ntz);
                    float[][] dxi=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,(float)scale);
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
    private static String stage(int scale,float[] trans,float[][][][] mem,float[][][][] shape) {
        float[] start=compute(trans,mem,shape,scale,true);float e0=start[2],energy=e0,lambda=1.0f;
        float[] finalTransform=trans;boolean accepted=false;int rejected=0;
        StringBuilder trials=new StringBuilder("[");
        while(true) {
            float[] candidate=new float[6];
            for(int i=0;i<6;i++)candidate[i]=trans[i]+lambda/start[9+i]*start[3+i];
            energy=compute(candidate,mem,shape,scale,false)[2];
            if(trials.length()>1)trials.append(',');
            trials.append("{\"transform\":").append(vector(candidate)).append(",\"energy\":").append(jsonNumber(energy)).append('}');
            if(energy>e0){finalTransform=candidate;lambda*=FACTOR;accepted=true;break;}
            lambda/=FACTOR;rejected++;if(rejected>MAX_REJECTIONS)break;
        }
        float diff=(energy-e0)/energy;
        String json="{\"scale\":"+scale+",\"input_dims\":"+vector(new float[]{mem.length,mem[0].length,mem[0][0].length})
            +",\"template_dims\":"+vector(new float[]{shape[1].length,shape[1][0].length,shape[1][0][0].length})
            +",\"initial_transform\":"+vector(trans)+",\"trials\":"+trials.append(']')
            +",\"transform\":"+vector(finalTransform)+",\"energy\":"+Float.toString(energy)
            +",\"lambda\":"+Float.toString(lambda)+",\"relative_change\":"+Float.toString(diff)
            +",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}';
        return json;
    }
    private static String runPyramid(int levels) {
        float[][][][] mem=fullMemberships(), sourceShape=fullShapes(), coarseShape=new float[CLASSES][][][];
        float[] transform={0,0,0,0,0,0};StringBuilder history=new StringBuilder("[");boolean first=true;
        for(int power=levels-1;power>0;power--) {
            int scale=1<<power;
            float[][][][] halfMem=ImageFunctions.subsample(mem,N,N,N,CLASSES,scale);
            for(int k=0;k<CLASSES;k++)coarseShape[k]=ImageFunctions.subsample(sourceShape[k],N,N,N,scale);
            String row=stage(scale,transform,halfMem,coarseShape);
            if(!first)history.append(',');first=false;history.append(row);
            // Parse the last transform from the local stage using its own calculation.
            float[] start=compute(transform,halfMem,coarseShape,scale,true);
            float e0=start[2],energy=e0,lambda=1.0f;float[] updated=transform;int rejected=0;boolean accepted=false;
            while(true){float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=transform[i]+lambda/start[9+i]*start[3+i];energy=compute(trial,halfMem,coarseShape,scale,false)[2];if(energy>e0){updated=trial;accepted=true;break;}lambda/=FACTOR;rejected++;if(rejected>MAX_REJECTIONS)break;}
            if(accepted)transform=updated;
        }
        return "{\"levels\":"+levels+",\"full_dimensions\":"+vector(new float[]{N,N,N})+",\"initial_transform\":[0,0,0,0,0,0],\"stages\":"+history.append(']')+",\"final_transform\":"+vector(transform)+'}';
    }
    public static void main(String[] args) {System.out.println("{\"cases\":["+runPyramid(1)+","+runPyramid(2)+","+runPyramid(3)+","+runPyramid(4)+"]}");}
}
