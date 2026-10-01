import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/** P18cye pinned source excerpt for the alternate approx=4 computeApproxMemberships() constructor path. */
public final class ToadsDefaultApproximateMembershipPassHarness {
    private static final String SEGMENTATION_SHA="63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String ALGORITHM_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String WRAPPER_SHA="d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115";
    private static final String DESCRIPTOR_SHA="d43efa5a002b0c1a838449ba7ca1d3bf35fab5e3b575474a71aadfe475dcc257";
    private static final String NUMERICS_SHA="13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb";
    private static final String POWERTABLE_SHA="5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final int CLASSES=11,APPROX=4;
    private static final float ZERO=1e-10f,INF=1e10f;

    private static final class Volume {
        final int nx,ny,nz,nc; final float[] f;
        Volume(Path path)throws IOException{
            byte[] raw;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(path));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);raw=out.toByteArray();}
            ByteBuffer h=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);if(h.getInt(0)!=348||h.getShort(70)!=16)throw new IOException("expected float32 NIfTI");
            nx=h.getShort(42);ny=h.getShort(44);nz=h.getShort(46);nc=h.getShort(48)>1?h.getShort(48):1;int start=(int)h.getFloat(108);f=new float[nx*ny*nz*nc];for(int i=0;i<f.length;i++)f[i]=h.getFloat(start+4*i);
        }
        int i(int x,int y,int z,int c){return x+nx*(y+ny*(z+nz*c));}
        float v(int x,int y,int z,int c){return f[i(x,y,z,c)];}
    }
    private static final class Table {
        final double q,step=(double)1.0e-6f; final double[] cache;
        Table(float exponent){q=exponent;int n=(int)Math.ceil(4.0/step)+1;cache=new double[n];for(int i=0;i<n;i++)cache[i]=Math.pow(i*step,q);}
        double lookup(float p){if(p<0||p>4)return Math.pow(p,q);int at=(int)(((double)p)/step+0.5);return cache[at];}
    }
    private static int[] bits(float[] a){int[] out=new int[a.length];for(int i=0;i<a.length;i++)out[i]=Float.floatToIntBits(a[i]);return out;}
    private static int[] ids(byte[] a){int[] out=new int[a.length];for(int i=0;i<a.length;i++)out[i]=a[i]&255;return out;}
    private static String json(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static float[] row(Path p)throws IOException{String[] t=Files.readString(p).trim().split("\\s+");float[] a=new float[t.length];for(int i=0;i<t.length;i++)a[i]=Float.parseFloat(t[i]);return a;}
    private static int flat(Volume v,int x,int y,int z,int c){return v.i(x,y,z,c);}
    private static final class Result {float[] imageVar,imageVarNum,imageVarDen,rawImageVar,imageStd,similarity,mems,den;byte[] best;float objective;}

    private static Result run(Volume image,Volume field,Volume mem,Volume bestV,Volume shape,Volume direct,
            float[] centers,float[] std,float[] classWeight,float[] exponent,float[] modality,float[] smoothPrior,float[] shapePrior,
            float priorScale,float priorFactor,float priorCoefficient,float smoothing,boolean useField){
        int nx=image.nx,ny=image.ny,nz=image.nz,n=nx*ny*nz;Result r=new Result();r.mems=new float[n*APPROX];r.best=new byte[n*APPROX];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){int vi=flat(image,x,y,z,0);for(int a=0;a<APPROX;a++){r.mems[vi*APPROX+a]=mem.v(x,y,z,a);r.best[vi*APPROX+a]=(byte)bestV.v(x,y,z,a);}}
        r.imageVar=new float[1];r.imageVarNum=new float[1];r.imageVarDen=new float[1];r.rawImageVar=new float[1];r.imageStd=new float[1];float num=0,den=0;
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int a=0;a<APPROX;a++){
            int k=(int)bestV.v(x,y,z,a);if(k==0)continue;float value=image.v(x,y,z,0);if(useField)value=field.v(x,y,z,0)*value;
            float diff=value-centers[k];diff=diff*diff;float m=mem.v(x,y,z,a);float weight=m*m+1e-30f;num+=weight*diff*diff;den+=weight;
        }
        float rawVar=num>0?den/num:Float.POSITIVE_INFINITY;r.imageVarNum[0]=num;r.imageVarDen[0]=den;r.rawImageVar[0]=rawVar;r.imageStd[0]=(float)Math.sqrt(rawVar);float sum=rawVar;r.imageVar[0]=rawVar/sum;
        // Selected one-channel normalization: rawVar / sum(rawVar) is exactly one.
        float[] sim=new float[CLASSES*CLASSES];for(int k=0;k<CLASSES;k++)for(int m=0;m<CLASSES;m++)if(m!=k){float d=0;float delta=centers[k]-centers[m];d+=r.imageVar[0]*delta*delta/(std[k]+std[m])/(std[k]+std[m]);sim[k*CLASSES+m]=1.0f/(1.0f+(1-priorScale)/priorScale*d);}r.similarity=sim;
        float objective=0; r.den=new float[n];Table[] tables=new Table[CLASSES];for(int k=1;k<CLASSES;k++)tables[k]=new Table(exponent[k]);
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int vi=flat(image,x,y,z,0);float[] scores=new float[CLASSES],energy=new float[CLASSES];float voxelShape[]=new float[CLASSES];for(int k=0;k<CLASSES;k++)voxelShape[k]=priorCoefficient>0?shape.v(x,y,z,k):0;
            float localDen=0;
            for(int k=0;k<CLASSES;k++){
                boolean isDirect=direct.v(x,y,z,k)!=0;if(!isDirect&&voxelShape[k]==0){scores[k]=0;continue;}
                float img=image.f[vi],likelihood;
                if(k==0){likelihood=modality[0]*classWeight[k]*r.imageVar[0]*(img==0?0:1);}
                else {float corrected=useField?field.f[vi]*img:img;float residual=corrected-centers[k];float distance=residual*residual;float weight=modality[0]*classWeight[k]*r.imageVar[0];likelihood=(float)(weight*tables[k].lookup(distance));}
                energy[k]=likelihood;float cost=likelihood;
                if(smoothing>0){float ngb=0;int[][] off={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};for(int[] o:off){int qi=flat(image,x+o[0],y+o[1],z+o[2],0);for(int a=0;a<APPROX;a++)if((r.best[qi*APPROX+a]&255)!=k){float m=r.mems[qi*APPROX+a];ngb+=m*m;}}cost+=smoothPrior[k]*smoothing*ngb;energy[k]+=0.5f*smoothPrior[k]*smoothing*ngb;}
                if(priorCoefficient>0){float inv=0;for(int m=0;m<CLASSES;m++)if(m!=k)inv+=shapePrior[m]*sim[k*CLASSES+m]*voxelShape[m]*voxelShape[m];cost+=priorCoefficient*inv;energy[k]+=0.5f*priorCoefficient*inv;}
                float value=cost>ZERO?1.0f/cost:INF;if(!isDirect)value=Math.min(1.0f,2.0f*priorFactor*voxelShape[k])*value;scores[k]=value;localDen+=value;
            }
            r.den[vi]=localDen;float[] scratch=scores.clone();for(int a=0;a<APPROX;a++){int winner=0;for(int k=1;k<CLASSES;k++)if(scratch[k]>scratch[winner])winner=k;r.best[vi*APPROX+a]=(byte)winner;r.mems[vi*APPROX+a]=scratch[winner];scratch[winner]=-INF;}
            if(localDen>0)for(int a=0;a<APPROX;a++){int k=r.best[vi*APPROX+a]&255;r.mems[vi*APPROX+a]/=localDen;objective+=energy[k]*r.mems[vi*APPROX+a]*r.mems[vi*APPROX+a];}
        }
        r.objective=objective;return r;
    }
    private static int at(Volume v,int x,int y,int z,int c){return v.i(x,y,z,c);}
    private static float[] orderedDen(Volume image,Result r){float[] out=new float[r.den.length];int n=0;for(int x=0;x<image.nx;x++)for(int y=0;y<image.ny;y++)for(int z=0;z<image.nz;z++)out[n++]=r.den[flat(image,x,y,z,0)];return out;}
    private static int[] edgeBest(Result r,Volume v,int x,int y,int z){int[] a=new int[APPROX];int i=at(v,x,y,z,0)*APPROX;for(int k=0;k<APPROX;k++)a[k]=r.best[i+k]&255;return a;}
    public static void main(String[] args)throws Exception{
        Volume image=new Volume(Path.of(args[0])),field=new Volume(Path.of(args[1])),mem=new Volume(Path.of(args[2])),best=new Volume(Path.of(args[3])),shape=new Volume(Path.of(args[4])),direct=new Volume(Path.of(args[5]));
        float[] center=row(Path.of(args[6])),std=row(Path.of(args[7])),vw=row(Path.of(args[8])),exp=row(Path.of(args[9])),mod=row(Path.of(args[10])),sp=row(Path.of(args[11])),ap=row(Path.of(args[12])),params=row(Path.of(args[13]));
        Result f=run(image,field,mem,best,shape,direct,center,std,vw,exp,mod,sp,ap,params[0],params[1],params[2],params[3],true);Result raw=run(image,field,mem,best,shape,direct,center,std,vw,exp,mod,sp,ap,params[0],params[1],params[2],params[3],false);
        int[] ids=edgeBest(f,image,9,9,9),zero=edgeBest(f,image,8,8,8),near=edgeBest(f,image,7,7,7),empty=edgeBest(f,image,10,10,10);
        System.out.print("{\"validation_profile\":\"alternate approx=4 constructor path; MedicAlgorithmToads08 shorter constructor forces useApprox=false\",\"segmentation_sha256\":\""+SEGMENTATION_SHA+"\",\"algorithm_sha256\":\""+ALGORITHM_SHA+"\",\"wrapper_sha256\":\""+WRAPPER_SHA+"\",\"descriptor_sha256\":\""+DESCRIPTOR_SHA+"\",\"numerics_sha256\":\""+NUMERICS_SHA+"\",\"powertable_sha256\":\""+POWERTABLE_SHA+"\",\"shape\":["+image.nx+","+image.ny+","+image.nz+",1],\"classes\":11,\"approx\":4,\"field_image_variance_bits\":"+json(bits(f.imageVar))+",\"field_image_variance_num_bits\":"+json(bits(f.imageVarNum))+",\"field_image_variance_den_bits\":"+json(bits(f.imageVarDen))+",\"field_image_variance_raw_bits\":"+json(bits(f.rawImageVar))+",\"field_image_std_bits\":"+json(bits(f.imageStd))+",\"raw_image_variance_bits\":"+json(bits(raw.imageVar))+",\"raw_image_variance_num_bits\":"+json(bits(raw.imageVarNum))+",\"raw_image_variance_den_bits\":"+json(bits(raw.imageVarDen))+",\"raw_image_variance_raw_bits\":"+json(bits(raw.rawImageVar))+",\"raw_image_std_bits\":"+json(bits(raw.imageStd))+",\"field_similarity_bits\":"+json(bits(f.similarity))+",\"field_membership_bits\":"+json(bits(f.mems))+",\"field_best\":"+json(ids(f.best))+",\"field_denominator_bits\":"+json(bits(orderedDen(image,f)))+",\"field_objective_bits\":"+Float.floatToIntBits(f.objective)+",\"raw_membership_bits\":"+json(bits(raw.mems))+",\"raw_objective_bits\":"+Float.floatToIntBits(raw.objective)+",\"tie_best\":"+json(ids)+",\"zero_energy_best\":"+json(zero)+",\"near_zero_best\":"+json(near)+",\"no_eligible_best\":"+json(empty)+",\"tie_voxel\":[9,9,9],\"zero_energy_voxel\":[8,8,8],\"near_zero_voxel\":[7,7,7],\"no_eligible_voxel\":[10,10,10]}");
    }
}
