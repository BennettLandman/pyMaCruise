import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/** Selected LesionToadAtlas.registerAllShapes() excerpt for P39zzzzt. */
public final class LesionRegistrationHandoffReferenceHarness {
    private static final int NX=6, NY=5, NZ=4, CLASSES=4;
    private static final float[] START={0.03f,-0.02f,0.01f,0.1f,-0.1f,0.05f};
    private static final float[] R={0.9f,1.2f,1.5f};
    private static final int SUBSAMPLE=2, IT_MAX=10, IT_SUPP=10;
    private static final float MIN_EDIFF=1e-6f, MIN_LAMBDA=0.001f, FACTOR=1.5f;
    private static final float[][][][] SHAPES=makeShapes();
    private static final ParametricTransform MODEL=new ParametricTransform(
        "rigid",NX/2.0f,NY/2.0f,NZ/2.0f,R[0],R[1],R[2],NX,NY,NZ,
        NX/2.0f,NY/2.0f,NZ/2.0f,R[0],R[1],R[2],NX,NY,NZ);

    private static float[][][][] makeShapes() {
        float[][][][] a=new float[CLASSES][NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            int f=(7*x+5*y+3*z+2*x*y+4*y*z)%11;
            a[0][x][y][z]=0.9f;
            a[1][x][y][z]=f<6?0.9f:0.02f;
            if(x==0)a[2][x][y][z]=0.95f;
            if(x==NX-1)a[2][x][y][z]=0.8f;
            if(y==0)a[2][x][y][z]=0.7f;
            if(y==NY-1)a[2][x][y][z]=0.6f;
        }
        return a;
    }
    private static float[][][][] loadMemberships(String file) throws Exception {
        float[][][][] mem=new float[NX][NY][NZ][CLASSES];
        List<Float> values=new ArrayList<>();
        try(BufferedReader r=new BufferedReader(new FileReader(file))) {
            String line;while((line=r.readLine())!=null)if(!line.trim().isEmpty())values.add(Float.parseFloat(line.trim()));
        }
        if(values.size()!=NX*NY*NZ*CLASSES)throw new IllegalArgumentException("membership count");
        int q=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<CLASSES;k++)mem[x][y][z][k]=values.get(q++);
        return mem;
    }
    private static float[] coefficients(float[][][][] mem,float[] trans,boolean derivative) {
        float[] h=new float[6],g=new float[6];float cost=0.0f,norm=0.0f;
        float[][] rot=ParametricTransform.computeRotation(trans);
        RotationMatrix rm=new RotationMatrix();rm.setParameters(trans[0],trans[1],trans[2]);
        float[][] da=rm.derivatives(1.0f,0.0f,0.0f),db=rm.derivatives(0.0f,1.0f,0.0f),dc=rm.derivatives(0.0f,0.0f,1.0f);
        for(int x=0;x<NX;x+=SUBSAMPLE)for(int y=0;y<NY;y+=SUBSAMPLE)for(int z=0;z<NZ;z+=SUBSAMPLE)for(int k=1;k<CLASSES;k++) {
            float[] xi=MODEL.imageToTemplate(x,y,z,trans,rot,1.0f);
            float p=ImageFunctions.linearInterpolation(SHAPES[k],0.0f,xi[0],xi[1],xi[2],NX,NY,NZ);
            if(p>0.0f) {
                float m=mem[x][y][z][k];
                if(derivative) {
                    float dx=ImageFunctions.linearInterpolationXderivative(SHAPES[k],xi[0],xi[1],xi[2],NX,NY,NZ);
                    float dy=ImageFunctions.linearInterpolationYderivative(SHAPES[k],xi[0],xi[1],xi[2],NX,NY,NZ);
                    float dz=ImageFunctions.linearInterpolationZderivative(SHAPES[k],xi[0],xi[1],xi[2],NX,NY,NZ);
                    float[][] dxi=MODEL.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,1.0f);
                    for(int i=0;i<6;i++) {
                        float dp=dx*dxi[0][i]+dy*dxi[1][i]+dz*dxi[2][i];
                        h[i]+=2.0f*m*m*dp*dp;
                        g[i]+=2.0f*m*m*p*dp;
                    }
                    cost+=m*m*p*p;
                } else cost+=m*m*p*p;
                norm+=p*p;
            }
        }
        float[] out=new float[14];out[0]=cost/norm;out[1]=norm;
        for(int i=0;i<6;i++){out[2+i]=g[i]/norm;out[8+i]=h[i]/norm;}
        return out;
    }
    private static String v(float[] a) {StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}return s.append(']').toString();}
    private static String matrix(float[][] a) {StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(v(a[i]));}return s.append(']').toString();}
    private static String bounds(float[] trans) {
        float[][] rot=ParametricTransform.computeRotation(trans);
        StringBuilder out=new StringBuilder("[");
        for(int k=0;k<CLASSES;k++) {
            if(k>0)out.append(',');int[] lo={NX,NY,NZ},hi={0,0,0};int count=0;
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
                float[] p=MODEL.imageToTemplate(x,y,z,trans,rot,1.0f);
                if(ImageFunctions.linearInterpolation(SHAPES[k],0.0f,p[0],p[1],p[2],NX,NY,NZ)>0.0f) {
                    count++;if(x<lo[0])lo[0]=x;if(y<lo[1])lo[1]=y;if(z<lo[2])lo[2]=z;
                    if(x>hi[0])hi[0]=x;if(y>hi[1])hi[1]=y;if(z>hi[2])hi[2]=z;
                }
            }
            out.append("{\"lower\":[").append(lo[0]).append(',').append(lo[1]).append(',').append(lo[2]).append("],\"upper\":[").append(hi[0]).append(',').append(hi[1]).append(',').append(hi[2]).append("],\"support_count\":").append(count).append('}');
        }
        return out.append(']').toString();
    }
    private static String optimize(String label,String file) throws Exception {
        float[][][][] mem=loadMemberships(file);float[] trans=START.clone();float lambda=1.0f,diff=1.0f;
        float energy=coefficients(mem,trans,false)[0];StringBuilder history=new StringBuilder("[");int n=0;
        while(n<IT_MAX && diff>MIN_EDIFF && lambda>MIN_LAMBDA) {
            if(n>0)history.append(',');float[] c=coefficients(mem,trans,true);float start=c[0],trialEnergy=start;int rejects=0;boolean accepted=false;
            StringBuilder trials=new StringBuilder("[");float[] acceptedTransform=trans;
            while(true) {
                float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=trans[i]+lambda/c[8+i]*c[2+i];
                trialEnergy=coefficients(mem,trial,false)[0];if(trials.length()>1)trials.append(',');trials.append("{\"transform\":").append(v(trial)).append(",\"energy\":").append(Float.toString(trialEnergy)).append('}');
                if(trialEnergy>start){acceptedTransform=trial;lambda*=FACTOR;accepted=true;break;}
                lambda/=FACTOR;rejects++;if(rejects>IT_SUPP)break;
            }
            float relative=(trialEnergy-start)/trialEnergy;diff=relative;
            if(accepted)trans=acceptedTransform;energy=trialEnergy;
            history.append("{\"trials\":").append(trials.append(']')).append(",\"accepted\":").append(accepted).append(",\"rejected_trials\":").append(rejects).append(",\"transform\":").append(v(trans)).append(",\"energy\":").append(Float.toString(energy)).append(",\"relative_change\":").append(Float.toString(relative)).append(",\"lambda\":").append(Float.toString(lambda)).append('}');n++;
        }
        String stop=n>=IT_MAX?"iteration_limit":(diff<=MIN_EDIFF?"min_energy_change":"lambda_floor");
        return "{\"name\":\""+label+"\",\"iterations\":"+history.append(']')+",\"transform\":"+v(trans)+",\"energy\":"+Float.toString(energy)+",\"lambda\":"+Float.toString(lambda)+",\"stop_reason\":\""+stop+"\",\"rotation\":"+matrix(ParametricTransform.computeRotation(trans))+",\"bounds\":"+bounds(trans)+"}";
    }
    public static void main(String[] args) throws Exception {
        System.out.println("{\"cases\":["+optimize("no_field",args[0])+","+optimize("with_field",args[1])+"]}");
    }
}
