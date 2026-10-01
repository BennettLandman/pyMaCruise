import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import edu.jhmi.rad.medic.structures.RotationMatrix;
import edu.jhmi.rad.medic.utilities.Numerics;

/**
 * Selected-source excerpt for initialApproxMemberships plus
 * registerShapesApproxPyramid. It uses vendored ParametricTransform,
 * RotationMatrix, ImageFunctions, and Numerics classes; no external checkout
 * is needed by the test suite.
 */
public final class ToadsDefaultApproximateRegistrationReferenceHarness {
    private static final String WRAPPER_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String ATLAS_SHA="e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2";
    private static final String DESCRIPTOR_SHA="cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939";
    private static final int APPROX=3, LEVELS=3, MAX_ITER=3, SUBSAMPLE=2;
    private static final float MIN_EDIFF=1e-6f, MIN_LAMBDA=1e-3f, FACTOR=1.5f;
    private static final float[] IMAGE_RES={.9f,1.2f,1.45f}, TEMPLATE_RES={.9f,1.2f,1.45f};
    private static final float[] INITIAL={.004f,-.003f,.002f,.04f,-.03f,.02f};

    private static final class Volume {
        final int nx,ny,nz,nc; final float[][][][] data;
        Volume(Path file) throws IOException {
            byte[] bytes;
            try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(file));ByteArrayOutputStream out=new ByteArrayOutputStream()){
                in.transferTo(out);bytes=out.toByteArray();
            }
            ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
            int dims=h.getShort(40),datatype=h.getShort(70);float offset=h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dims>3?h.getShort(48):1;
            if(h.getInt(0)!=348||datatype!=16||nx<1||ny<1||nz<1||nc<1)throw new IOException("expected float32 NIfTI-1");
            data=new float[nx][ny][nz][nc];int start=(int)offset;
            for(int c=0;c<nc;c++)for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
                int i=x+nx*(y+ny*(z+nz*c));data[x][y][z][c]=h.getFloat(start+4*i);
            }
        }
    }
    private static final class Level {
        int scale,nx,ny,nz,tx,ty,tz;float[][][][] mem,shape;byte[][][][] best;
        float[] initial,finalTransform;float initialEnergy,energy,lambda,difference;
        String stopReason;
        List<String> iterations=new ArrayList<>();
    }
    private static int bits(float x){return Float.floatToIntBits(x);}
    private static String intArray(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static String byteArray(byte[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(byte[][][] x:a)for(byte[][] y:x)for(byte[] z:y)for(byte v:z){if(!first)s.append(',');first=false;s.append(v&255);}return s.append(']').toString();}
    private static String floatArray(float[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){if(!first)s.append(',');first=false;s.append(bits(v));}return s.append(']').toString();}
    private static String vector(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String floats(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}return s.append(']').toString();}
    private static List<String> classTypes(Path descriptor)throws IOException{
        List<String> lines=Files.readAllLines(descriptor);int head=-1;
        for(int i=0;i<lines.size();i++)if(lines.get(i).startsWith("Structures:")){head=i;break;}
        if(head<0)throw new IOException("descriptor lacks Structures header");
        int count=Integer.parseInt(lines.get(head).split(":",2)[1].trim());List<String> types=new ArrayList<>();
        for(int k=0;k<count;k++){
            String[] fields=lines.get(head+1+k).trim().split("\\s+");
            String tag=fields[fields.length-1];
            if(tag.equals("mask"))types.add("MASK");else if(tag.equals("optimized_csf"))types.add("OPTIMIZED_CSF");
            else if(tag.equals("optimized_gm"))types.add("OPTIMIZED_GM");
            else if(tag.equals("optimized_wm"))types.add("OPTIMIZED_WM");
            else if(tag.equals("lesion"))types.add("LESION");else types.add("OBJECT");
        }
        return types;
    }
    private static float[][] centroids(int channels,int classes){
        float[][] out=new float[channels][classes];for(int c=0;c<channels;c++)for(int k=0;k<classes;k++)out[c][k]=.19f+.021f*k+.013f*c;return out;
    }
    private static void initialize(Volume image,float[][] center,float[][][][] mem,byte[][][][] best){
        int classes=center[0].length,channels=image.nc;
        for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
            float[] num=new float[classes],grouping=new float[classes];float den=0.0f;
            for(int k=0;k<classes;k++){
                num[k]=0.0f;
                for(int c=0;c<channels;c++){
                    float delta=image.data[x][y][z][c]-center[c][k];num[k]+=delta*delta;
                }
                if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
                den+=num[k];
                grouping[k]=0.0f;for(int m=0;m<classes;m++)if(center[0][m]==center[0][k])grouping[k]++;
            }
            if(den>0.0f)for(int k=0;k<classes;k++)num[k]=grouping[k]*num[k]/den;
            Numerics.bestIndex(best[x][y][z],mem[x][y][z],num,APPROX);
        }
    }
    private static Level pool(float[][][][] mem,byte[][][][] best,int classes,int factor){
        Level q=new Level();q.scale=factor;q.nx=mem.length/factor;q.ny=mem[0].length/factor;q.nz=mem[0][0].length/factor;
        q.mem=new float[q.nx][q.ny][q.nz][APPROX];q.best=new byte[q.nx][q.ny][q.nz][APPROX];float div=factor*factor*factor;
        for(int x=0;x<q.nx;x++)for(int y=0;y<q.ny;y++)for(int z=0;z<q.nz;z++){
            float[] val=new float[classes];
            for(int i=0;i<factor;i++)for(int j=0;j<factor;j++)for(int l=0;l<factor;l++)for(int a=0;a<APPROX;a++){
                int k=best[x*factor+i][y*factor+j][z*factor+l][a];val[k]+=mem[x*factor+i][y*factor+j][z*factor+l][a]/div;
            }
            Numerics.bestIndex(q.best[x][y][z],q.mem[x][y][z],val,APPROX);
        }
        return q;
    }
    private static float[][][][] downsampleShapes(Volume atlas,int factor){
        float[][][][] out=new float[atlas.nc][atlas.nx/factor][atlas.ny/factor][atlas.nz/factor];
        for(int k=0;k<atlas.nc;k++){
            float[][][] one=new float[atlas.nx][atlas.ny][atlas.nz];
            for(int x=0;x<atlas.nx;x++)for(int y=0;y<atlas.ny;y++)for(int z=0;z<atlas.nz;z++)one[x][y][z]=atlas.data[x][y][z][k];
            out[k]=ImageFunctions.subsample(one,atlas.nx,atlas.ny,atlas.nz,factor);
        }
        return out;
    }
    private static float[] mapped(ParametricTransform model,int x,int y,int z,float[] trans,float[][] rot,float scale){return model.imageToTemplate(x,y,z,trans,rot,scale);}
    private static float measure(Level l,float[] trans,ParametricTransform model,boolean derivatives,float[] gradOut,float[] hessOut){
        RotationMatrix r=model.computeRotationMatrix(trans);float[][] rot=r.getMatrix();
        float[][] da=r.derivatives(1,0,0),db=r.derivatives(0,1,0),dc=r.derivatives(0,0,1);
        float cost=0.0f,norm=0.0f;float[] grad=new float[6],hess=new float[6];
        for(int x=0;x<l.nx;x+=SUBSAMPLE)for(int y=0;y<l.ny;y+=SUBSAMPLE)for(int z=0;z<l.nz;z+=SUBSAMPLE)for(int a=0;a<APPROX;a++){
            int k=l.best[x][y][z][a];if(k<=0)continue;
            float[] pxyz=mapped(model,x,y,z,trans,rot,l.scale);
            float prior=ImageFunctions.linearClosestInterpolation(l.shape[k],pxyz[0],pxyz[1],pxyz[2],l.tx,l.ty,l.tz);
            if(prior<=0.0f)continue;
            float m=l.mem[x][y][z][a],m2=m*m,p2=prior*prior;
            cost+=m2*p2;norm+=p2;
            if(derivatives){
                float px=ImageFunctions.linearInterpolationXderivative(l.shape[k],pxyz[0],pxyz[1],pxyz[2],l.tx,l.ty,l.tz);
                float py=ImageFunctions.linearInterpolationYderivative(l.shape[k],pxyz[0],pxyz[1],pxyz[2],l.tx,l.ty,l.tz);
                float pz=ImageFunctions.linearInterpolationZderivative(l.shape[k],pxyz[0],pxyz[1],pxyz[2],l.tx,l.ty,l.tz);
                float[][] jac=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,l.scale);float[] dp=new float[6];
                for(int i=0;i<6;i++)dp[i]=px*jac[0][i]+py*jac[1][i]+pz*jac[2][i];
                for(int i=0;i<6;i++){grad[i]+=2.0f*m2*prior*dp[i];hess[i]+=2.0f*m2*dp[i]*dp[i];}
            }
        }
        if(derivatives){for(int i=0;i<6;i++){gradOut[i]=grad[i]/norm;hessOut[i]=hess[i]/norm;}}
        return cost/norm;
    }
    private static Level runLevel(Level l,float[] trans,int[] fullImage,int[] fullTemplate,float[] imageRes,float[] templateRes){
        ParametricTransform model=new ParametricTransform("rigid",fullImage[0]/2.0f,fullImage[1]/2.0f,fullImage[2]/2.0f,
            imageRes[0],imageRes[1],imageRes[2],fullImage[0],fullImage[1],fullImage[2],
            fullTemplate[0]/2.0f,fullTemplate[1]/2.0f,fullTemplate[2]/2.0f,
            templateRes[0],templateRes[1],templateRes[2],fullTemplate[0],fullTemplate[1],fullTemplate[2]);
        l.initial=trans.clone();float lambda=1.0f,diff=1.0f,energy=measure(l,trans,model,false,null,null);l.initialEnergy=energy;
        for(int n=0;n<MAX_ITER&&diff>MIN_EDIFF&&lambda>MIN_LAMBDA;n++){
            float[] grad=new float[6],hess=new float[6];float e0=measure(l,trans,model,true,grad,hess),e=e0;
            StringBuilder trials=new StringBuilder("[");int rejected=0;boolean accepted=false;float[] candidate=trans.clone();
            while(true){
                float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=trans[i]+lambda/Math.max(Numerics.ZERO,hess[i])*grad[i];
                e=measure(l,trial,model,false,null,null);if(trials.length()>1)trials.append(',');
                trials.append("{\"lambda\":").append(bits(lambda)).append(",\"transform\":").append(vector(trial)).append(",\"energy\":").append(bits(e)).append('}');
                if(e>e0){candidate=trial;lambda=lambda*FACTOR;accepted=true;break;}
                lambda=lambda/FACTOR;rejected++;if(rejected>10)break;
            }
            trans=candidate;diff=(e-e0)/e;energy=e;
            l.iterations.add("{\"gradient\":"+vector(grad)+",\"hessian\":"+vector(hess)+",\"trials\":"+trials.append(']')+
                ",\"transform\":"+vector(trans)+",\"energy\":"+bits(energy)+",\"lambda\":"+bits(lambda)+",\"relative_change\":"+bits(diff)+",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}');
        }
        if(diff<=MIN_EDIFF)l.stopReason="energy_tolerance";
        else if(lambda<=MIN_LAMBDA)l.stopReason="lambda_floor";
        else l.stopReason="iteration_limit";
        l.finalTransform=trans.clone();l.energy=energy;l.lambda=lambda;l.difference=diff;return l;
    }
    private static String levelJson(Level l){
        return "{\"scale\":"+l.scale+",\"dims\":["+l.nx+","+l.ny+","+l.nz+"],\"template_dims\":["+l.tx+","+l.ty+","+l.tz+"],"+
            "\"memberships_bits\":"+floatArray(l.mem)+",\"best_indices\":"+byteArray(l.best)+",\"initial_transform_bits\":"+vector(l.initial)+
            ",\"initial_energy_bits\":"+bits(l.initialEnergy)+",\"iterations\":["+String.join(",",l.iterations)+"],\"transform_bits\":"+vector(l.finalTransform)+
            ",\"energy_bits\":"+bits(l.energy)+",\"lambda_bits\":"+bits(l.lambda)+",\"relative_change_bits\":"+bits(l.difference)+
            ",\"stop_reason\":\""+l.stopReason+"\"}";
    }
    public static void main(String[] args)throws Exception{
        List<String> types=classTypes(Path.of(args[2]));
        if(!types.equals(Arrays.asList("MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_GM","OPTIMIZED_WM","OPTIMIZED_WM","OPTIMIZED_WM")))
            throw new IllegalArgumentException("unexpected default TOADS atlas descriptor order");
        Volume image=new Volume(Path.of(args[0])),atlas=new Volume(Path.of(args[1]));
        if(image.nc!=2||atlas.nc!=11)throw new IllegalArgumentException("fixture geometry/classes mismatch");
        int classes=atlas.nc;float[][][][] initMem=new float[image.nx][image.ny][image.nz][APPROX];byte[][][][] initBest=new byte[image.nx][image.ny][image.nz][APPROX];
        float[][] centers=centroids(image.nc,classes);centers[0][2]=centers[0][1];
        initialize(image,centers,initMem,initBest);
        List<String> stageJson=new ArrayList<>();float[] trans=INITIAL.clone();int[] fi={image.nx,image.ny,image.nz},ft={atlas.nx,atlas.ny,atlas.nz};
        for(int scale=1<<(LEVELS-1);scale>=2;scale/=2){
            Level l=pool(initMem,initBest,classes,scale);l.shape=downsampleShapes(atlas,scale);l.tx=l.shape[0].length;l.ty=l.shape[0][0].length;l.tz=l.shape[0][0][0].length;
            l=runLevel(l,trans,fi,ft,IMAGE_RES,TEMPLATE_RES);trans=l.finalTransform;stageJson.add(levelJson(l));
        }
        StringBuilder out=new StringBuilder("{\"wrapper_source_sha256\":\"").append(WRAPPER_SHA).append("\",\"preprocess_source_sha256\":\"").append(PREPROCESS_SHA)
            .append("\",\"atlas_source_sha256\":\"").append(ATLAS_SHA).append("\",\"descriptor_sha256\":\"").append(DESCRIPTOR_SHA).append("\",\"approx\":").append(APPROX).append(",\"levels\":").append(LEVELS)
            .append(",\"subsample\":").append(SUBSAMPLE).append(",\"initial_transform_bits\":").append(vector(INITIAL))
            .append(",\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"]")
            .append(",\"call_order\":[\"initialMemberships\",\"setImageInfo\",\"initShapeRegistration\",\"setTransform(preprocess)\",\"updateShapeRegistration\",\"registerShapesPyramid\"]")
            .append(",\"pyramid_dispatch\":\"registerShapesApproxPyramid\",\"scale_order\":[4,2],\"coarse_subsample\":2")
            .append(",\"initial_dims\":[").append(image.nx).append(',').append(image.ny).append(',').append(image.nz)
            .append("],\"initial_memberships_bits\":").append(floatArray(initMem)).append(",\"initial_best_indices\":").append(byteArray(initBest))
            .append(",\"registered_class_indices\":[1,2,3,4,5,6,7,8,9,10],\"stages\":[").append(String.join(",",stageJson))
            .append("],\"final_transform_bits\":").append(vector(trans)).append('}');
        System.out.println(out);
    }
}
