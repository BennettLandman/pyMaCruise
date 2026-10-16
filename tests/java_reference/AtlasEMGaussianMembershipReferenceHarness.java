/*
 * One-channel/no-field Atlas-mode computeMemberships excerpt from the pinned
 * AtlasEMSegmentation source. Atlas priors/covariances are fixed test inputs;
 * this harness omits registration and classification output generation.
 */
public final class AtlasEMGaussianMembershipReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final float PI=3.14159265358979f, Log2PI=(float)Math.log(PI), SMOOTH=0.625f;
    static final int NX=4,NY=4,NZ=4,K=3;
    static final float[][][][] image=new float[NX][NY][NZ][1];
    static final float[][][][] mem=new float[NX][NY][NZ][K];
    static final float[][][][] prior=new float[NX][NY][NZ][K];
    static final float[][] mean={{0.15f,0.47f,0.83f}};
    static final float[][][] covariance={{{0.12f}},{{0.2f}},{{0.31f}}};
    static final boolean[][][] mask=new boolean[NX][NY][NZ];
    static void init(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            mask[x][y][z]=x>0&&x<NX-1&&y>0&&y<NY-1&&z>0&&z<NZ-1;
            image[x][y][z][0]=(float)(x*7+y*5+z*3+1)/20.0f;
            for(int k=0;k<K;k++){
                mem[x][y][z][k]=(float)(1+x+y*2+z*3+k)/13.0f;
                prior[x][y][z][k]=mask[x][y][z]?(float)(1+x*2+y+z*3+k)/20.0f:0.0f;
            }
        }
    }
    static float[] update(){
        float[] old=new float[K],energy=new float[K]; float Energy=0.0f,distance=0.0f;
        float[][] inverse=new float[K][1]; float[] det=new float[K];
        for(int k=0;k<K;k++){det[k]=Math.abs(covariance[k][0][0]);inverse[k][0]=1.0f/covariance[k][0][0];}
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]){
            float den=1e-30f;
            for(int k=0;k<K;k++){
                old[k]=mem[x][y][z][k];
                if(prior[x][y][z][k]!=0)
                    energy[k]=(float)(-(1/2.0f)*Log2PI-0.5*Math.log(det[k])+Math.log(prior[x][y][z][k]));
                else energy[k]=0;
                float delta=image[x][y][z][0]-mean[0][k];
                float exponent=0.0f;
                float temp=inverse[k][0]*delta;
                exponent+=temp*delta;
                exponent/=2.0f;
                float ngb=0.0f;
                int[][] n={{x+1,y,z},{x-1,y,z},{x,y+1,z},{x,y-1,z},{x,y,z+1},{x,y,z-1}};
                for(int[] q:n)for(int m=0;m<K;m++)if(m!=k)ngb+=mem[q[0]][q[1]][q[2]][m];
                exponent+=SMOOTH*ngb;
                energy[k]-=exponent;
                mem[x][y][z][k]=(float)(Math.exp(-exponent)*prior[x][y][z][k]/Math.sqrt(det[k]));
                den+=mem[x][y][z][k];
            }
            for(int k=0;k<K;k++)mem[x][y][z][k]=mem[x][y][z][k]/den;
            for(int k=0;k<K;k++){
                Energy+=(energy[k]+Math.log(mem[x][y][z][k])*old[k]);
                float change=Math.abs(mem[x][y][z][k]-old[k]);if(change>distance)distance=change;
            }
        }else for(int k=0;k<K;k++)mem[x][y][z][k]=0;
        return new float[]{Energy,distance};
    }
    static String vector(float[] a){String s="[";for(int i=0;i<a.length;i++){if(i>0)s+=",";s+=Float.floatToRawIntBits(a[i]);}return s+"]";}
    static String tensor(){String s="[";boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){if(!first)s+=",";first=false;s+=Float.floatToRawIntBits(mem[x][y][z][k]);}return s+"]";}
    public static void main(String[] args){init();float[] metrics=update();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"memberships\":"+tensor()+",\"metrics\":"+vector(metrics)+"}");}
}
