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
 * Selected-source excerpt for the atlas-aware 13-class approximate-membership
 * overload only. The wrapper uses this overload after segmentation state
 * exists, then updates and registers shapes at that later stage; this harness
 * does not claim it is the wrapper's initial pyramid call.
 */
public final class ToadsLesionApproximateRegistrationReferenceHarness {
    private static final String WRAPPER_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String LESION_WRAPPER_SHA="556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String ATLAS_SHA="e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2";
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
            String tag=lines.get(head+1+k).split("\\t")[2];
            if(tag.equals("mask"))types.add("MASK");else if(tag.equals("optimized_csf"))types.add("OPTIMIZED_CSF");
            else if(tag.equals("lesion"))types.add("LESION");else types.add("OBJECT");
        }
        return types;
    }
    private static float[][] centroids(int channels,int classes){
        float[][] out=new float[channels][classes];for(int c=0;c<channels;c++)for(int k=0;k<classes;k++)out[c][k]=.19f+.021f*k+.013f*c;return out;
    }
    private static Object[] initializeVoxel(float[] image,float[][] center,float[] shape){
        int classes=shape.length;float[] score=new float[classes];float den=0.0f;
        for(int k=0;k<classes;k++){
            score[k]=0.0f;
            for(int c=0;c<image.length;c++){float delta=image[c]-center[c][k];score[k]+=delta*delta;}
            for(int m=0;m<classes;m++)if(m!=k)score[k]+=shape[m]*shape[m];
            if(score[k]>1e-20f)score[k]=1.0f/score[k];else score[k]=1e20f;
            den+=score[k];
        }
        if(den>0.0f)for(int k=0;k<classes;k++)score[k]/=den;
        byte[] best=new byte[APPROX];float[] selected=new float[APPROX];
        Numerics.bestIndex(best,selected,score,APPROX,Numerics.bestIndex(shape));
        return new Object[]{best,selected};
    }
    private static String byteVector(byte[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]&255);}return s.append(']').toString();}
    private static void initialize(Volume image,Volume prior,float[][] center,float[][][][] mem,byte[][][][] best){
        int classes=prior.nc,channels=image.nc;
        for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
            float[] num=new float[classes];float den=0.0f;
            for(int k=0;k<classes;k++){
                num[k]=0.0f;
                for(int c=0;c<channels;c++){
                    float delta=image.data[x][y][z][c]-center[c][k];num[k]+=delta*delta;
                }
                for(int m=0;m<classes;m++)if(m!=k){float value=prior.data[x][y][z][m];num[k]+=value*value;}
                if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
                den+=num[k];
            }
            if(den>0.0f)for(int k=0;k<classes;k++)num[k]/=den;
            float[] localShape=new float[classes];for(int k=0;k<classes;k++)localShape[k]=prior.data[x][y][z][k];
            Numerics.bestIndex(best[x][y][z],mem[x][y][z],num,APPROX,Numerics.bestIndex(localShape));
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
        l.finalTransform=trans.clone();l.energy=energy;l.lambda=lambda;l.difference=diff;return l;
    }
    private static String levelJson(Level l){
        return "{\"scale\":"+l.scale+",\"dims\":["+l.nx+","+l.ny+","+l.nz+"],\"template_dims\":["+l.tx+","+l.ty+","+l.tz+"],"+
            "\"memberships_bits\":"+floatArray(l.mem)+",\"best_indices\":"+byteArray(l.best)+",\"initial_transform_bits\":"+vector(l.initial)+
            ",\"initial_energy_bits\":"+bits(l.initialEnergy)+",\"iterations\":["+String.join(",",l.iterations)+"],\"transform_bits\":"+vector(l.finalTransform)+
            ",\"energy_bits\":"+bits(l.energy)+",\"lambda_bits\":"+bits(l.lambda)+",\"relative_change_bits\":"+bits(l.difference)+'}';
    }
    public static void main(String[] args)throws Exception{
        List<String> types=classTypes(Path.of(args[3]));
        if(!types.equals(Arrays.asList("MASK","OPTIMIZED_CSF","OPTIMIZED_CSF","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","OBJECT","LESION","OBJECT","OBJECT")))
            throw new IllegalArgumentException("unexpected lesion-bearing atlas descriptor order");
        Volume image=new Volume(Path.of(args[0])),subjectShapes=new Volume(Path.of(args[1])),atlas=new Volume(Path.of(args[2]));
        if(image.nx!=subjectShapes.nx||image.ny!=subjectShapes.ny||image.nz!=subjectShapes.nz||image.nc!=2||subjectShapes.nc!=13||atlas.nc!=13)throw new IllegalArgumentException("fixture geometry/classes mismatch");
        int classes=atlas.nc;float[][][][] initMem=new float[image.nx][image.ny][image.nz][APPROX];byte[][][][] initBest=new byte[image.nx][image.ny][image.nz][APPROX];
        initialize(image,subjectShapes,centroids(image.nc,classes),initMem,initBest);
        StringBuilder out=new StringBuilder("{\"wrapper_source_sha256\":\"").append(WRAPPER_SHA).append("\",\"preprocess_source_sha256\":\"").append(PREPROCESS_SHA)
            .append("\",\"atlas_source_sha256\":\"").append(ATLAS_SHA).append("\",\"approx\":").append(APPROX)
            .append(",\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OBJECT\",\"OBJECT\",\"OBJECT\",\"OBJECT\",\"OBJECT\",\"OBJECT\",\"OBJECT\",\"LESION\",\"OBJECT\",\"OBJECT\"]")
            .append(",\"validation_scope\":\"method-only atlas-aware 13-class approximate-membership overload\"")
            .append(",\"selected_overload\":\"preprocess.initialMemberships(segmentation.getMemberships(),segmentation.getBest(),atlas,approx,centroids,classes)\"")
            .append(",\"later_wrapper_context\":[\"segmentation state exists\",\"initialMemberships(...,atlas,...)\",\"updateShapeRegistration\",\"registerShapes\"]")
            .append(",\"initial_dims\":[").append(image.nx).append(',').append(image.ny).append(',').append(image.nz)
            .append("],\"initial_memberships_bits\":").append(floatArray(initMem)).append(",\"initial_best_indices\":").append(byteArray(initBest)).append('}');
        out.setLength(out.length()-1);
        out.append(",\"lesion_wrapper_source_sha256\":\"").append(LESION_WRAPPER_SHA).append('\"')
            .append(",\"selected_wrapper\":\"AlgorithmMultipleToads\"")
            .append(",\"source_call_order\":[\"new MultipleToadSegmentation(...)\",\"segmentation.initCentroids(centroids)\",\"segmentation.setLesionCentroid(...)\",\"preprocess.initialMemberships(segmentation.getMemberships(),segmentation.getBest(),atlas,approx,centroids,classes)\",\"atlas.updateShapeRegistration(...)\",\"atlas.registerShapes() when alignShapes\"]")
            .append(",\"reachability_classification\":{\"selected\":\"atlas-aware overload is reached from AlgorithmMultipleToads after segmentation initialization\",\"initial_pyramid\":\"earlier AlgorithmMultipleToads call is the no-atlas overload\",\"dedicated_lesion_wrapper\":\"AlgorithmLesionToads uses the no-atlas overload for its initial registration pyramid\",\"scope\":\"overload method only; no registration result claimed\"}");
        float[] tieImage={0.0f};float[][] tieCentroids=new float[1][13];tieCentroids[0][5]=1.0f;float[] tieShape=new float[13];tieShape[5]=1.0f;
        Object[] tie=initializeVoxel(tieImage,tieCentroids,tieShape);
        out.append(",\"shape_seeded_membership_tie\":{\"indices\":").append(byteVector((byte[])tie[0]))
            .append(",\"membership_bits\":").append(vector((float[])tie[1])).append('}');
        Object[] zero=initializeVoxel(new float[]{0.0f},new float[1][13],new float[13]);
        out.append(",\"zero_cost_tie\":{\"indices\":").append(byteVector((byte[])zero[0]))
            .append(",\"membership_bits\":").append(vector((float[])zero[1])).append('}').append('}');
        System.out.println(out);
    }
}
