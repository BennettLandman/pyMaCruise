import edu.jhmi.rad.medic.utilities.PowerTable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** Selected Java formula harness on explicit next-loop state arrays (not a full TOADS run). */
public final class ToadsP18czpActualNextLoopMembershipReferenceHarness {
    static final int N=14,C=11,NC=2,V=N*N*N;
    static final int[] CATEGORY={-1,0,0,1,1,1,1,1,2,2,2};
    static final float SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f, PRIOR_SCALE=.05f;
    static final float[] MODALITY={.73f,.46f};
    static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    static final float[][] VARIANCE_WEIGHT={
        {.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},
        {1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}
    };
    static final float[][] EXPONENT=new float[NC][C];
    final byte[] labels=new byte[V],closest=new byte[V*4];
    final float[] images=new float[V*NC],fields=new float[V*NC],mem=new float[V*C],shapes=new float[V*C],old=new float[NC*C];
    final float[] centers=new float[NC*C],deviations=new float[NC*C],imageVariance=new float[NC],imageStddev=new float[NC];
    final float[] similarity=new float[C*C],energy=new float[V*C],denominator=new float[V],updated=new float[V*C];
    float objective;
    int direct,neighbor,shapeOnly;

    int idx(int x,int y,int z){return (x*N+y)*N+z;}
    void read(Path path)throws Exception {
        byte[] raw=Files.readAllBytes(path); ByteBuffer b=ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        int expected=V+V*NC*4*2+V*C*4*2+NC*C*4+V*4;
        if(raw.length!=expected)throw new IllegalArgumentException("unexpected explicit-state byte count: "+raw.length+" != "+expected);
        b.get(labels);
        for(int i=0;i<images.length;i++)images[i]=b.getFloat();
        for(int i=0;i<fields.length;i++)fields[i]=b.getFloat();
        for(int i=0;i<mem.length;i++)mem[i]=b.getFloat();
        for(int i=0;i<old.length;i++)old[i]=b.getFloat();
        for(int i=0;i<shapes.length;i++)shapes[i]=b.getFloat();
        b.get(closest);
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)EXPONENT[c][k]=CATEGORY[k]<0?1.0f:(CATEGORY[k]==0?.8f:1.0f);
    }
    void statistics(){
        PowerTable[] tables={new PowerTable(0,4,.000001f,.8f),new PowerTable(0,4,.000001f,1.0f),new PowerTable(0,4,.000001f,1.0f)};
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++){
            float num=0,den=0;
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
                int p=idx(x,y,z);float corrected=fields[p*NC+c]*images[p*NC+c],m=mem[p*C+k],w=m*m;
                if(k>0){float residual=corrected-old[c*C+k],d=residual*residual,f=1;
                    if(d>1e-10f)f=(float)(tables[CATEGORY[k]].lookup((double)d,(double)EXPONENT[c][k])/d);
                    w=f*m;w*=m;
                }
                num+=w*corrected;den+=w;
            }
            centers[c*C+k]=den>0?.9f*num/den+.1f*old[c*C+k]:old[c*C+k];
        }
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++){
            float num=0,den=0;
            for(int p=0;p<V;p++)if(labels[p]==k){float m=mem[p*C+k],w=m*m+1e-30f;
                float residual=fields[p*NC+c]*images[p*NC+c]-centers[c*C+k];num+=w*residual*residual;den+=w;}
            deviations[c*C+k]=den>0?(float)Math.sqrt(num/den):0;
        }
        float total=0;
        for(int c=0;c<NC;c++){
            float num=0,den=0;
            for(int k=1;k<C;k++)for(int p=0;p<V;p++){
                float m=mem[p*C+k],w=m*m+1e-30f,res=fields[p*NC+c]*images[p*NC+c]-centers[c*C+k],d=res*res;
                num+=w*d*d;den+=w;
            }
            imageVariance[c]=num>0?den/num:1e10f;total+=imageVariance[c];
        }
        for(int c=0;c<NC;c++){imageStddev[c]=(float)Math.sqrt(imageVariance[c]);imageVariance[c]/=total;}
        float scale=(1-PRIOR_SCALE)/PRIOR_SCALE;
        for(int k=0;k<C;k++)for(int m=0;m<C;m++)if(k!=m){float distance=0;
            for(int c=0;c<NC;c++){float d=centers[c*C+k]-centers[c*C+m],s=deviations[c*C+k]+deviations[c*C+m];float t=imageVariance[c]*d;t*=d;t/=s;t/=s;distance+=t;}
            similarity[k*C+m]=1/(1+scale*distance);
        }
    }
    void memberships(){
        PowerTable[][] tables=new PowerTable[3][NC];for(int t=0;t<3;t++)for(int c=0;c<NC;c++)tables[t][c]=new PowerTable(0,4,.000001f,t==0?.8f:1.0f);
        System.arraycopy(mem,0,updated,0,mem.length);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<C;k++){
            int p=idx(x,y,z);float e=0;
            if(CATEGORY[k]<0){for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=imageVariance[c];if(images[p*NC+c]!=0)e+=w;}}
            else for(int c=0;c<NC;c++){float residual=fields[p*NC+c]*images[p*NC+c]-centers[c*C+k],d=residual*residual;
                double term=tables[CATEGORY[k]][c].lookup((double)d,(double)EXPONENT[c][k]);float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=imageVariance[c];e=(float)(e+w*term);}
            energy[p*C+k]=e;
        }
        int[][] off={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            int p=idx(x,y,z),base=p*4;float den=0;float[] obj=new float[C];
            for(int k=0;k<C;k++){
                boolean directClass=closest[base]==k||closest[base+1]==k,near=closest[base+2]==k||closest[base+3]==k;
                if(!directClass&&!near&&shapes[p*C+k]==0){energy[p*C+k]=0;updated[p*C+k]=0;continue;}
                if(directClass)direct++;else if(near)neighbor++;else shapeOnly++;
                float num=energy[p*C+k],o=num,ns=0;
                for(int[] d:off)for(int m=0;m<C;m++)if(m!=k){float v=updated[idx(x+d[0],y+d[1],z+d[2])*C+m];ns+=v*v;}
                num+=SMOOTH_PRIOR[k]*SMOOTH*ns;o+=.5f*SMOOTH_PRIOR[k]*SMOOTH*ns;
                float inv=0;for(int m=0;m<C;m++)if(m!=k){float v=ATLAS_PRIOR[m]*similarity[k*C+m];v*=shapes[p*C+m];v*=shapes[p*C+m];inv+=v;}
                num+=COEFFICIENT*inv;o+=.5f*COEFFICIENT*inv;
                float f=shapes[p*C+k];if(directClass)f=1;else if(near)f=Math.max(.25f,Math.min(2*FACTOR*f,1));else f=Math.min(1,2*FACTOR*f);
                float val=f*(num>1e-10f?1/num:1e10f);updated[p*C+k]=val;den+=val;obj[k]=o;
            }
            denominator[p]=den;if(den>0)for(int k=0;k<C;k++){updated[p*C+k]/=den;float v=obj[k]*updated[p*C+k];objective+=v*updated[p*C+k];}
        }
    }
    static String hash(float[] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float v:a){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}StringBuilder s=new StringBuilder();for(byte v:d.digest())s.append(String.format("%02x",v));return s.toString();}
    static String bits(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(a[i]));}return s.append(']').toString();}
    public static void main(String[] args)throws Exception{
        if(args.length!=1)throw new IllegalArgumentException("usage: harness explicit-state.bin");
        ToadsP18czpActualNextLoopMembershipReferenceHarness h=new ToadsP18czpActualNextLoopMembershipReferenceHarness();h.read(Path.of(args[0]));h.statistics();h.memberships();
        System.out.print("{\"reference_scope\":\"selected Java formula harness on explicitly supplied next-loop state; fixed two-channel parameters and 11-class descriptor; no outer TOADS iteration\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"wrapper_source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"power_table_source_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"centroids_bits\":"+bits(h.centers)+",\"deviations_bits\":"+bits(h.deviations)+",\"image_variance_bits\":"+bits(h.imageVariance)+",\"image_stddev_bits\":"+bits(h.imageStddev)+",\"similarity_bits\":"+bits(h.similarity)+",\"membership_sha256\":\""+hash(h.updated)+"\",\"denominator_sha256\":\""+hash(h.denominator)+"\",\"class_energy_sha256\":\""+hash(h.energy)+"\",\"objective_bits\":"+Float.floatToIntBits(h.objective)+",\"shape_branch_counts\":["+h.direct+","+h.neighbor+","+h.shapeOnly+"]}");
    }
}
