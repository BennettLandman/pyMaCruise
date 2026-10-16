/* generateClassification excerpt from AtlasEMSegmentation.java. */
public final class AtlasEMClassificationReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int NX=3,NY=2,NZ=2,K=3;
    static final float[][][][] mem=new float[NX][NY][NZ][K];
    static final byte[] labels={(byte)4,(byte)17,(byte)91};
    static byte[][][] classify(){
        byte[][][] out=new byte[NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            float max=0;int best=-1;
            for(int k=0;k<K;k++)if(mem[x][y][z][k]>max){best=k;max=mem[x][y][z][k];}
            if(best>-1)out[x][y][z]=labels[best];else out[x][y][z]=0;
        }
        return out;
    }
    static String tensor(byte[][][] a){String s="[";boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){if(!first)s+=",";first=false;s+=(int)a[x][y][z];}return s+"]";}
    public static void main(String[] args){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++)mem[x][y][z][k]=(x+y+z+k)%5==0?0.0f:(float)(1+x*2+y+z+k)/9.0f;
        mem[0][0][0][0]=0.5f;mem[0][0][0][1]=0.5f;mem[0][0][0][2]=0.1f;
        mem[1][1][1][0]=0.0f;mem[1][1][1][1]=0.0f;mem[1][1][1][2]=0.0f;
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"labels\":"+tensor(classify())+"}");
    }
}
