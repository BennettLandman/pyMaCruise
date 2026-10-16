/*
 * Fixed-prior no-field T1 AtlasEM wrapper loop excerpt. This tests the
 * source-order parameter/membership updates with synthetic transformed priors.
 */
public final class AtlasEMIterationReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int NX=5,NY=5,NZ=5,K=3,MAX=4;
    static final float THRESHOLD=0.01f,SMOOTH=0.125f;
    static final float[][][][] image=new float[NX][NY][NZ][1];
    static final float[][][][] mem=new float[NX][NY][NZ][K];
    static final float[][][][] prior=new float[NX][NY][NZ][K];
    static final boolean[][][] mask=new boolean[NX][NY][NZ];
    static final float[][] centroid={{0.2f,0.5f,0.8f}};
    static final float[][][] cov={{{0.2f}},{{0.3f}},{{0.4f}}};
    static final float[] energies=new float[MAX],changes=new float[MAX];
    static int iterations=0;static boolean converged=false;
    static void init(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            mask[x][y][z]=x>0&&x<NX-1&&y>0&&y<NY-1&&z>0&&z<NZ-1;
            image[x][y][z][0]=(float)(x*11+y*7+z*3+1)/80.0f;
            for(int k=0;k<K;k++){
                mem[x][y][z][k]=(float)(1+x+y*2+z*3+k)/23.0f;
                prior[x][y][z][k]=mask[x][y][z]?(float)(1+x*2+y+z*3+k)/20.0f:0.0f;
            }
        }
    }
    static void run(){
        for(int it=0;it<MAX;it++){
            for(int k=0;k<K;k++){
                float num=0,den=0;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]){num+=mem[x][y][z][k]*image[x][y][z][0];den+=mem[x][y][z][k];}
                if(den!=0)centroid[0][k]=num/den;
            }
            for(int k=0;k<K;k++){
                float num=0,den=0;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]){
                    float d=image[x][y][z][0]-centroid[0][k];num+=d*d*mem[x][y][z][k];den+=mem[x][y][z][k];
                }
                if(den!=0)cov[k][0][0]=num/den;
            }
            float[] old=new float[K];float distance=0,energy=0;
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]){
                float den=1e-30f;float[] loglike=new float[K];
                for(int k=0;k<K;k++){
                    old[k]=mem[x][y][z][k];float p=prior[x][y][z][k];
                    loglike[k]=p!=0?(float)(-(1/2.0f)*((float)Math.log(PI))-0.5*Math.log(Math.abs(cov[k][0][0]))+Math.log(p)):0.0f;
                    float delta=image[x][y][z][0]-centroid[0][k];float exp=0;
                    float temp=(1.0f/cov[k][0][0])*delta;exp+=temp*delta;exp/=2.0f;
                    if(SMOOTH>0){float ngb=0;int[][] n={{x+1,y,z},{x-1,y,z},{x,y+1,z},{x,y-1,z},{x,y,z+1},{x,y,z-1}};
                        for(int[] q:n)for(int m=0;m<K;m++)if(m!=k)ngb+=mem[q[0]][q[1]][q[2]][m];exp+=SMOOTH*ngb;}
                    loglike[k]-=exp;
                    mem[x][y][z][k]=(float)(Math.exp(-exp)*p/Math.sqrt(Math.abs(cov[k][0][0])));den+=mem[x][y][z][k];
                }
                for(int k=0;k<K;k++)mem[x][y][z][k]/=den;
                for(int k=0;k<K;k++){energy+=loglike[k]+Math.log(mem[x][y][z][k])*old[k];float diff=Math.abs(mem[x][y][z][k]-old[k]);if(diff>distance)distance=diff;}
            }else for(int k=0;k<K;k++)mem[x][y][z][k]=0;
            energies[it]=energy;changes[it]=distance;iterations=it+1;
            if(distance<THRESHOLD){converged=true;break;}
        }
    }
    static final float PI=3.14159265358979f;
    static String vector(float[] a,int n){String s="[";for(int i=0;i<n;i++){if(i>0)s+=",";s+=Float.floatToRawIntBits(a[i]);}return s+"]";}
    static String tensor(){String s="[";boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){if(!first)s+=",";first=false;s+=Float.floatToRawIntBits(mem[x][y][z][k]);}return s+"]";}
    static String matrix(){String s="[";for(int k=0;k<K;k++){if(k>0)s+=",";s+="["+Float.floatToRawIntBits(centroid[0][k])+"]";}return s+"]";}
    static String covariances(){String s="[";for(int k=0;k<K;k++){if(k>0)s+=",";s+="["+Float.floatToRawIntBits(cov[k][0][0])+"]";}return s+"]";}
    public static void main(String[] args){init();run();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"iterations\":"+iterations+",\"converged\":"+converged+",\"energies\":"+vector(energies,iterations)+",\"changes\":"+vector(changes,iterations)+",\"centroids\":"+matrix()+",\"covariances\":"+covariances()+",\"memberships\":"+tensor()+"}");}
}
