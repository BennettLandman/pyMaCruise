import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;

/** Selected-source formula reference for one-channel exact TOADS post-thinning stages. */
public final class ToadsPostThinningExactIterationReferenceHarness {
    static final int C=11;
    static final float[] MODALITY={.73f};
    static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    static final float[] VARIANCE_WEIGHT={.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f};
    static final float SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f, PRIOR_SCALE=.05f;

    static int bits(float value){return Float.floatToRawIntBits(value);}
    static String ints(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    static String hash(float[] values)throws Exception{
        MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        for(float v:values){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}
        StringBuilder s=new StringBuilder();for(byte v:d.digest())s.append(String.format("%02x",v));return s.toString();
    }
    static int at(int x,int y,int z,int ny,int nz){return (x*ny+y)*nz+z;}
    static String intHash(int[] values)throws Exception{float[] f=new float[values.length];for(int i=0;i<f.length;i++)f[i]=Float.intBitsToFloat(values[i]);return hash(f);}

    public static void main(String[] args)throws Exception{
        int nx=Integer.parseInt(args[1]),ny=Integer.parseInt(args[2]),nz=Integer.parseInt(args[3]);
        int volume=nx*ny*nz;
        byte[] labels=new byte[volume],closest=new byte[volume*4];
        float[] image=new float[volume],memberships=new float[volume*C],field=new float[volume],old=new float[C];
        try(DataInputStream in=new DataInputStream(new FileInputStream(args[0]))){
            for(int i=0;i<volume;i++)labels[i]=in.readByte();
            for(int i=0;i<volume;i++)image[i]=in.readFloat();
            for(int i=0;i<volume*C;i++)memberships[i]=in.readFloat();
            for(int i=0;i<volume;i++)field[i]=in.readFloat();
            for(int k=0;k<C;k++)old[k]=in.readFloat();
            in.readFully(closest);
        }
        float[] center=old.clone(),std=new float[C];
        PowerTable csfCentroid=new PowerTable(0.0f,4.0f,.000001f,.8f);
        PowerTable tissueCentroid=new PowerTable(0.0f,4.0f,.000001f,1.0f);
        for(int k=0;k<C;k++){
            float num=0.0f,den=0.0f;
            for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
                int v=at(x,y,z,ny,nz),m=v*C+k;float corrected=field[v]*image[v],mem=memberships[m],weight;
                if(k==0)weight=mem*mem;
                else{
                    float residual=corrected-old[k],distance=residual*residual,factor=1.0f;
                    if(distance>1e-10f){PowerTable t=k<=2?csfCentroid:tissueCentroid;float exponent=k<=2?.8f:1.0f;factor=(float)(t.lookup((double)distance,(double)exponent)/(double)distance);}
                    weight=factor*mem;weight*=mem;
                }
                num+=weight*corrected;den+=weight;
            }
            if(den>0.0f)center[k]=.9f*num/den+.1f*old[k];
        }
        for(int k=0;k<C;k++){
            float num=0.0f,den=0.0f;
            for(int v=0;v<volume;v++)if(labels[v]==k){float mem=memberships[v*C+k],weight=mem*mem+1e-30f;
                float corrected=field[v]*image[v],diff=corrected-center[k];num+=weight*diff*diff;den+=weight;}
            if(den>0.0f)std[k]=(float)Math.sqrt(num/den);
        }
        float imageNum=0.0f,imageDen=0.0f;
        for(int k=1;k<C;k++)for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){
            int v=at(x,y,z,ny,nz);float mem=memberships[v*C+k],weight=mem*mem+1e-30f;
            float corrected=field[v]*image[v],residual=corrected-center[k],diff=residual*residual;
            float term=weight*diff;term*=diff;imageNum+=term;imageDen+=weight;
        }
        float imageVariance=imageNum>0.0f?imageDen/imageNum:1e10f;
        float imageStd=(float)Math.sqrt(imageVariance);
        imageVariance=imageVariance/imageVariance;
        float[][] similarity=new float[C][C];
        float scale=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;
        for(int k=0;k<C;k++)for(int m=0;m<C;m++)if(k!=m){float diff=center[k]-center[m],sum=std[k]+std[m];float term=imageVariance*diff;term*=diff;term/=sum;term/=sum;similarity[k][m]=1.0f/(1.0f+scale*term);}

        float[] energy=new float[volume*C],denominator=new float[volume];
        double[] powerAudit=new double[volume*C];
        PowerTable[] powers={new PowerTable(0.0f,4.0f,.000001f,.8f),new PowerTable(0.0f,4.0f,.000001f,1.0f),new PowerTable(0.0f,4.0f,.000001f,1.0f)};
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int v=at(x,y,z,ny,nz);
            for(int k=0;k<C;k++){
                if(k==0){if(image[v]!=0.0f)energy[v*C+k]=MODALITY[0]*VARIANCE_WEIGHT[k]*imageVariance;}
                else{float corrected=field[v]*image[v],residual=corrected-center[k],d=residual*residual;int category=k<=2?0:k<=7?1:2;
                    double exponent=category==0?(double).8f:1.0;
                    double term=powers[category].lookup((double)d,exponent);
                    powerAudit[v*C+k]=term;
                    float w=MODALITY[0]*VARIANCE_WEIGHT[k];w*=imageVariance;energy[v*C+k]=(float)(w*term);}
            }
        }
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        float objective=0.0f;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int v=at(x,y,z,ny,nz);float[] objectiveAt=new float[C];float den=0.0f;
            for(int k=0;k<C;k++){
                boolean direct=false,neighbor=false;for(int i=0;i<2;i++)if(closest[v*4+i]==k)direct=true;for(int i=2;i<4;i++)if(closest[v*4+i]==k)neighbor=true;
                if(!direct&&!neighbor){energy[v*C+k]=0.0f;memberships[v*C+k]=0.0f;continue;}
                float num=energy[v*C+k],obj=num,neighborSum=0.0f;
                for(int[] d:offsets)for(int m=0;m<C;m++)if(m!=k){float q=memberships[at(x+d[0],y+d[1],z+d[2],ny,nz)*C+m];neighborSum+=q*q;}
                float increment=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;num+=increment;obj+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                float inverse=0.0f; // zero transformed-shape inputs make inverse atlas prior zero.
                num+=COEFFICIENT*inverse;obj+=.5f*COEFFICIENT*inverse;
                float factor=direct?1.0f:.25f;
                float value=factor*(num>1e-10f?1.0f/num:1e10f);memberships[v*C+k]=value;den+=value;objectiveAt[k]=obj;
            }
            denominator[v]=den;
            if(den>0.0f)for(int k=0;k<C;k++){memberships[v*C+k]/=den;float contribution=objectiveAt[k]*memberships[v*C+k];objective+=contribution*memberships[v*C+k];}
        }
        int[] centerBits=new int[C],stdBits=new int[C],similarityBits=new int[C*C];
        for(int k=0;k<C;k++){centerBits[k]=bits(center[k]);stdBits[k]=bits(std[k]);for(int m=0;m<C;m++)similarityBits[k*C+m]=bits(similarity[k][m]);}
        try(DataOutputStream out=new DataOutputStream(new FileOutputStream(args[4]))){for(float value:energy)out.writeFloat(value);for(float value:memberships)out.writeFloat(value);for(float value:denominator)out.writeFloat(value);}
        try(DataOutputStream out=new DataOutputStream(new FileOutputStream(args[5]))){for(double value:powerAudit)out.writeDouble(value);}
        String harnessHash="formula_operations_from_selected_source_not_source_class_execution";
        System.out.print("{\"reference_kind\":\"source_formula_harness\",\"source_order\":[\"computeCorrectionField.degree1\",\"computeExactCentroids\",\"computeExactVariances\",\"computeImageVariance\",\"centroid similarity\",\"computeExactMemberships\"],\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"power_table_source_sha256\":\"5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678\",\"shape\":["+nx+","+ny+","+nz+"],\"centroid_bits\":"+ints(centerBits)+",\"standard_deviation_bits\":"+ints(stdBits)+",\"image_variance_bits\":"+bits(imageVariance)+",\"image_std_bits\":"+bits(imageStd)+",\"similarity_bits\":"+ints(similarityBits)+",\"membership_sha256\":\""+hash(memberships)+"\",\"denominator_sha256\":\""+hash(denominator)+"\",\"class_energy_sha256\":\""+hash(energy)+"\",\"objective_bits\":"+bits(objective)+"}");
    }
}
