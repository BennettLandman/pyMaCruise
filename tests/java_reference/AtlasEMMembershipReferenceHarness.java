/* Fixed-state AtlasEMSegmentation.computeMemberships excerpt: Atlas priors,
 * one T1 channel, no correction field. Source loop and neighbor order retained. */
public final class AtlasEMMembershipReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int N=4,K=3;
    static final float[][][][] image=new float[N][N][N][1], mem=new float[N][N][N][K], prior=new float[N][N][N][K];
    static final boolean[][][] mask=new boolean[N][N][N];
    static final float[] mean={0.19f,0.61f,0.84f},cov={0.08f,0.13f,0.21f};
    static final float smooth=10.0f/6.0f/K;
    static int bits(float x){return Float.floatToRawIntBits(x);}
    static void init(){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            mask[x][y][z]=x>0&&x<N-1&&y>0&&y<N-1&&z>0&&z<N-1;
            image[x][y][z][0]=(float)((x*19+y*7+z*3)%17)/16.0f;
            for(int k=0;k<K;k++){
                mem[x][y][z][k]=(1+x+y*2+z*3+k)/19.0f;
                prior[x][y][z][k]=(x+y+z+k)%7==0?0.0f:(1+x*2+y+z+k*3)/23.0f;
            }
        }
    }
    static Object[] update(){
        float distance=0.0f, Energy=0.0f;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)if(mask[x][y][z]){
            float[] old=mem[x][y][z].clone(), energy=new float[K]; float den=1e-30f;
            for(int k=0;k<K;k++){
                float p=prior[x][y][z][k];
                if(p!=0.0f)energy[k]=(float)(-0.5f*(float)Math.log(3.14159265358979f)-0.5*Math.log(Math.abs(cov[k]))+Math.log(p));
                else energy[k]=0.0f;
                if(p>0.0f){
                    float d=image[x][y][z][0]-mean[k];
                    float exponent=(1.0f/cov[k])*d*d; exponent/=2.0f;
                    float ngb=0.0f;
                    int[][] q={{x+1,y,z},{x-1,y,z},{x,y+1,z},{x,y-1,z},{x,y,z+1},{x,y,z-1}};
                    for(int[] v:q)for(int m=0;m<K;m++)if(m!=k)ngb+=mem[v[0]][v[1]][v[2]][m];
                    exponent+=smooth*ngb; energy[k]-=exponent;
                    mem[x][y][z][k]=(float)(Math.exp(-exponent)*p/Math.sqrt(Math.abs(cov[k])));
                    den+=mem[x][y][z][k];
                }else mem[x][y][z][k]=0.0f;
            }
            if(den>0.0f)for(int k=0;k<K;k++)mem[x][y][z][k]/=den;else for(int k=0;k<K;k++)mem[x][y][z][k]=0.0f;
            for(int k=0;k<K;k++){
                float change=Math.abs(mem[x][y][z][k]-old[k]);if(change>distance)distance=change;
                Energy+=(energy[k]+Math.log(mem[x][y][z][k])*old[k]);
            }
        }else for(int k=0;k<K;k++)mem[x][y][z][k]=0.0f;
        int[] label=new int[N*N*N];int p=0;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int best=0;for(int k=1;k<K;k++)if(mem[x][y][z][k]>mem[x][y][z][best])best=k;
            label[p++]=mem[x][y][z][best]>0.0f?best:-1;
        }
        float[] flat=new float[N*N*N*K];p=0;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<K;k++)flat[p++]=mem[x][y][z][k];
        return new Object[]{flat,label,Energy,distance};
    }
    static String ints(int[] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=a[i];}return s+"]";}
    static String bits(float[] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=Integer.toUnsignedString(Float.floatToRawIntBits(a[i]));}return s+"]";}
    public static void main(String[] args){init();Object[] r=update();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"memberships\":"+bits((float[])r[0])+",\"labels\":"+ints((int[])r[1])+
        ",\"energy\":"+Integer.toUnsignedString(bits((float)r[2]))+",\"distance\":"+Integer.toUnsignedString(bits((float)r[3]))+"}");}
}
