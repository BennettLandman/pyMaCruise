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
import edu.jhmi.rad.medic.utilities.Numerics;

/** Pinned-method reference for atlas-aware top-k membership and scale-one rigid registration. */
public final class ToadsDefaultSingleScaleRegistrationReferenceHarness {
    private static final String WRAPPER_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String PREPROCESS_SHA="31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7";
    private static final String ATLAS_SHA="e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2";
    private static final int APPROX=3, SUBSAMPLE=2, MAX_ITER=1, IT_SUPP=10;
    private static final float MIN_EDIFF=1e-6f, MIN_LAMBDA=1e-3f, FACTOR=1.5f;
    private static final float[] RES={.9f,1.2f,1.45f}, INITIAL={.031f,-.027f,.022f,.31f,-.27f,.22f};
    private static final class Volume {
        final int nx,ny,nz,nc; final float[][][][] data;
        Volume(Path file)throws IOException{
            byte[] bytes;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(file));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);bytes=out.toByteArray();}
            ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);int dims=h.getShort(40),datatype=h.getShort(70);float offset=h.getFloat(108);
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=dims>3?h.getShort(48):1;
            if(h.getInt(0)!=348||datatype!=16)throw new IOException("expected float32 NIfTI-1");
            data=new float[nx][ny][nz][nc];int start=(int)offset;
            for(int c=0;c<nc;c++)for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++)data[x][y][z][c]=h.getFloat(start+4*(x+nx*(y+ny*(z+nz*c))));
        }
    }
    private static int bits(float x){return Float.floatToIntBits(x);}
    private static String vector(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String floats(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.toString(a[i]));}return s.append(']').toString();}
    private static String flat(float[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){if(!first)s.append(',');first=false;s.append(bits(v));}return s.append(']').toString();}
    private static String flat(byte[][][][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(byte[][][] x:a)for(byte[][] y:x)for(byte[] z:y)for(byte v:z){if(!first)s.append(',');first=false;s.append(v&255);}return s.append(']').toString();}
    private static String matrix(float[][] a){StringBuilder s=new StringBuilder("[");boolean first=true;for(float[] row:a)for(float v:row){if(!first)s.append(',');first=false;s.append(bits(v));}return s.append(']').toString();}
    private static float[][] centers(int channels){float[][] c=new float[channels][11];for(int ch=0;ch<channels;ch++)for(int k=0;k<11;k++)c[ch][k]=.19f+.021f*k+.013f*ch;return c;}
    private static void initialize(Volume image,Volume prior,float[][] center,float[][][][] mem,byte[][][][] best){
        for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
            float[] num=new float[11],shape=new float[11];
            for(int k=0;k<11;k++){
                for(int ch=0;ch<image.nc;ch++){float d=image.data[x][y][z][ch]-center[ch][k];num[k]+=d*d;}
                for(int m=0;m<11;m++)if(m!=k){float p=prior.data[x][y][z][m];num[k]+=p*p;}
                if(num[k]>1e-20f)num[k]=1.0f/num[k];else num[k]=1e20f;
                shape[k]=prior.data[x][y][z][k];
            }
            float den=0;for(float v:num)den+=v;if(den>0)for(int k=0;k<11;k++)num[k]/=den;
            Numerics.bestIndex(best[x][y][z],mem[x][y][z],num,APPROX,Numerics.bestIndex(shape));
        }
    }
    private static String run(Volume image,Volume atlas,Volume prior){
        int[] fi={image.nx,image.ny,image.nz},ft={atlas.nx,atlas.ny,atlas.nz};
        float[][][][] mem=new float[image.nx][image.ny][image.nz][APPROX];byte[][][][] best=new byte[image.nx][image.ny][image.nz][APPROX];
        initialize(image,prior,centers(image.nc),mem,best);
        ParametricTransform model=new ParametricTransform("rigid",fi[0]/2f,fi[1]/2f,fi[2]/2f,RES[0],RES[1],RES[2],fi[0],fi[1],fi[2],ft[0]/2f,ft[1]/2f,ft[2]/2f,RES[0],RES[1],RES[2],ft[0],ft[1],ft[2]);
        float[] trans=INITIAL.clone();RotationMatrix rotation=model.computeRotationMatrix(trans);float[][] rot=rotation.getMatrix();float[][] da=rotation.derivatives(1,0,0),db=rotation.derivatives(0,1,0),dc=rotation.derivatives(0,0,1);
        List<String> iterations=new ArrayList<>();float lambda=1,diff=1,energy=0;
        for(int outer=0;outer<MAX_ITER&&diff>MIN_EDIFF&&lambda>MIN_LAMBDA;outer++){
            float cost=0,norm=0;float[] grad=new float[6],hess=new float[6];
            for(int x=0;x<image.nx;x+=SUBSAMPLE)for(int y=0;y<image.ny;y+=SUBSAMPLE)for(int z=0;z<image.nz;z+=SUBSAMPLE)for(int a=0;a<APPROX;a++){
                int k=best[x][y][z][a]&255;if(k<=0)continue;
                float[] xyz=model.imageToTemplate(x,y,z,trans,rot,1f);float[] shape=new float[atlas.nx*atlas.ny*atlas.nz];int ii=0;
                for(int sz=0;sz<atlas.nz;sz++)for(int sy=0;sy<atlas.ny;sy++)for(int sx=0;sx<atlas.nx;sx++)shape[ii++]=atlas.data[sx][sy][sz][k];
                float p=ImageFunctions.linearClosestInterpolation(shape,xyz[0],xyz[1],xyz[2],ft[0],ft[1],ft[2]);if(p<=0)continue;
                float m=mem[x][y][z][a],m2=m*m,p2=p*p;cost+=m2*p2;norm+=p2;
                float px=ImageFunctions.linearInterpolationXderivative(shape,xyz[0],xyz[1],xyz[2],ft[0],ft[1],ft[2]);
                float py=ImageFunctions.linearInterpolationYderivative(shape,xyz[0],xyz[1],xyz[2],ft[0],ft[1],ft[2]);
                float pz=ImageFunctions.linearInterpolationZderivative(shape,xyz[0],xyz[1],xyz[2],ft[0],ft[1],ft[2]);
                float[][] jac=model.imageToTemplateDerivatives(x,y,z,trans,rot,da,db,dc,1f);float[] dp=new float[6];
                for(int i=0;i<6;i++)dp[i]=px*jac[0][i]+py*jac[1][i]+pz*jac[2][i];
                for(int i=0;i<6;i++){grad[i]+=2f*m*m*p*dp[i];hess[i]+=2f*m*m*dp[i]*dp[i];}
            }
            float e0=cost/norm;for(int i=0;i<6;i++){grad[i]/=norm;hess[i]/=norm;}
            float[] coeffGrad=grad.clone(),coeffHess=hess.clone();StringBuilder trials=new StringBuilder("[");int rejected=0;boolean accepted=false;float e=e0;float[] acceptedTransform=trans.clone();
            while(true){float[] trial=new float[6];for(int i=0;i<6;i++)trial[i]=trans[i]+lambda/hess[i]*grad[i];
                float trialCost=0,trialNorm=0;float[][] trialRot=model.computeRotation(trial);
                for(int x=0;x<image.nx;x+=SUBSAMPLE)for(int y=0;y<image.ny;y+=SUBSAMPLE)for(int z=0;z<image.nz;z+=SUBSAMPLE)for(int a=0;a<APPROX;a++){
                    int k=best[x][y][z][a]&255;if(k<=0)continue;float[] xyz=model.imageToTemplate(x,y,z,trial,trialRot,1f);float[] sh=new float[atlas.nx*atlas.ny*atlas.nz];int ii=0;for(int sz=0;sz<atlas.nz;sz++)for(int sy=0;sy<atlas.ny;sy++)for(int sx=0;sx<atlas.nx;sx++)sh[ii++]=atlas.data[sx][sy][sz][k];float p=ImageFunctions.linearClosestInterpolation(sh,xyz[0],xyz[1],xyz[2],ft[0],ft[1],ft[2]);if(p>0){float m=mem[x][y][z][a];trialCost+=m*m*p*p;trialNorm+=p*p;}}
                e=trialCost/trialNorm;if(trials.length()>1)trials.append(',');trials.append("{\"lambda\":").append(bits(lambda)).append(",\"transform\":").append(vector(trial)).append(",\"energy\":").append(bits(e)).append('}');
                if(e>e0){acceptedTransform=trial;lambda*=FACTOR;accepted=true;break;}lambda/=FACTOR;rejected++;if(rejected>IT_SUPP)break;}
            diff=(e-e0)/e;trans=acceptedTransform;energy=e;
            iterations.add("{\"gradient\":"+vector(coeffGrad)+",\"hessian\":"+vector(coeffHess)+",\"trials\":"+trials.append(']')+",\"transform\":"+vector(trans)+",\"energy\":"+bits(energy)+",\"lambda\":"+bits(lambda)+",\"relative_change\":"+bits(diff)+",\"accepted\":"+accepted+",\"rejected_trials\":"+rejected+'}');
        }
        String stop=diff<=MIN_EDIFF?"energy_tolerance":lambda<=MIN_LAMBDA?"lambda_floor":"iteration_limit";
        float[][] finalRotation=model.computeRotation(trans);StringBuilder bounds=new StringBuilder("[");
        for(int k=0;k<11;k++){
            int[] lo={image.nx,image.ny,image.nz},hi={0,0,0};int count=0;float[] shape=new float[atlas.nx*atlas.ny*atlas.nz];int ii=0;
            for(int z=0;z<atlas.nz;z++)for(int y=0;y<atlas.ny;y++)for(int x=0;x<atlas.nx;x++)shape[ii++]=atlas.data[x][y][z][k];
            float[][] mat=new float[3][4];model.precomputeImageToTemplateMatrix(mat,trans,finalRotation,1f);
            for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++){
                float[] xyz={mat[0][0]*x+mat[0][1]*y+mat[0][2]*z+mat[0][3],mat[1][0]*x+mat[1][1]*y+mat[1][2]*z+mat[1][3],mat[2][0]*x+mat[2][1]*y+mat[2][2]*z+mat[2][3]};
                if(ImageFunctions.linearClosestInterpolation(shape,xyz[0],xyz[1],xyz[2],atlas.nx,atlas.ny,atlas.nz)>0){int[] pt={x,y,z};count++;for(int a=0;a<3;a++){lo[a]=Math.min(lo[a],pt[a]);hi[a]=Math.max(hi[a],pt[a]);}}
            }
            if(k>0)bounds.append(',');bounds.append('[').append(lo[0]).append(',').append(lo[1]).append(',').append(lo[2]).append(',').append(hi[0]).append(',').append(hi[1]).append(',').append(hi[2]).append(',').append(count).append(']');
        }
        bounds.append(']');
        return "{\"wrapper_source_sha256\":\""+WRAPPER_SHA+"\",\"preprocess_source_sha256\":\""+PREPROCESS_SHA+"\",\"atlas_source_sha256\":\""+ATLAS_SHA+"\",\"approx\":"+APPROX+",\"class_types\":[\"MASK\",\"OPTIMIZED_CSF\",\"OPTIMIZED_CSF\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_GM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\",\"OPTIMIZED_WM\"],\"call_order\":[\"initialApproxMemberships(atlas)\",\"updateShapeRegistration(rigid,iter=1,level=1)\",\"registerShapes()->registerAllShapes(single-scale)\"],\"initial_transform_bits\":"+vector(INITIAL)+",\"initial_memberships_bits\":"+flat(mem)+",\"best_indices\":"+flat(best)+",\"iterations\":["+String.join(",",iterations)+"],\"final_transform_bits\":"+vector(trans)+",\"rotation_bits\":"+matrix(finalRotation)+",\"bounds\":"+bounds+",\"final_energy_bits\":"+bits(energy)+",\"lambda_bits\":"+bits(lambda)+",\"stop_reason\":\""+stop+"\"}";
    }
    public static void main(String[] args)throws Exception{System.out.println(run(new Volume(Path.of(args[0])),new Volume(Path.of(args[1])),new Volume(Path.of(args[2]))));}
}
