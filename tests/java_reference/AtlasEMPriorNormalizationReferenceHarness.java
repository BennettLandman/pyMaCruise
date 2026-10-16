/* Selected computeAtlas atlas-prior normalization excerpt. */
public final class AtlasEMPriorNormalizationReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int NX=3,NY=2,NZ=2,K=3;
    static final float[][][][] prior=new float[K][NX][NY][NZ];
    static final boolean[][][] mask=new boolean[NX][NY][NZ];
    static void init(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            mask[x][y][z]=(x+y+z)%2==0;
            for(int k=0;k<K;k++)prior[k][x][y][z]=(float)(x*3+y*2+z+k)/7.0f;
        }
        prior[0][0][0][0]=0;prior[1][0][0][0]=0;prior[2][0][0][0]=0;
    }
    static void normalize(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(mask[x][y][z]){
            float norm=0;for(int k=0;k<K;k++)norm+=prior[k][x][y][z];
            if(norm!=0)for(int k=0;k<K;k++)prior[k][x][y][z]/=norm;
        }
    }
    static String tensor(){String s="[";boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++){if(!first)s+=",";first=false;s+=Float.floatToRawIntBits(prior[k][x][y][z]);}return s+"]";}
    public static void main(String[] args){init();normalize();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"priors\":"+tensor()+"}");}
}
