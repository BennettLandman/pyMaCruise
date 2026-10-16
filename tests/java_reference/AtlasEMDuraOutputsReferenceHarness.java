/* Fixed-input excerpt from AtlasEMSegmentation.generateDuraRemovalOutputs;
 * MIPAV FILL_HOLES border/flood behavior copied from AlgorithmMorphology25D. */
import java.util.ArrayDeque;
import java.util.Queue;
public final class AtlasEMDuraOutputsReferenceHarness {
    static final String SOURCE_SHA="b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a";
    static final int X=7,Y=7,Z=2,K=6;
    static final String[] names={"Background","SulcalCSF","CorticalGM","OtherA","Cerebellum-GM","OtherB"};
    static final float[][][][] mem=new float[X][Y][Z][K];
    static final int[][][] seg=new int[X][Y][Z];
    static final boolean[][][] mask=new boolean[X][Y][Z];
    static final float[][][][] out=new float[X][Y][Z][3];
    static void init(){
        for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++){
            mask[x][y][z]=(x+y+z)%5!=0;seg[x][y][z]=(x*3+y+z*2)%K;
            for(int k=0;k<K;k++)mem[x][y][z][k]=(1+x*3+y*5+z*7+k*11)/101.0f;
        }
    }
    static void create(){
        for(int k=0;k<K;k++)if(!names[k].equals("SulcalCSF")&&!names[k].equals("Sulcal-CSF")&&!names[k].equals("CorticalGM")&&!names[k].equals("Cerebrum-GM")&&!names[k].equals("CerebralGM")&&!names[k].equals("Cerebellum-GM")&&!names[k].equals("CerebellarGM")&&!names[k].equals("GM")&&!names[k].equals("Background"))
            for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)if(mask[x][y][z]){out[x][y][z][1]+=mem[x][y][z][k];if(seg[x][y][z]==k)out[x][y][z][2]=1.0f;}
        for(int k=0;k<K;k++)if(names[k].equals("CorticalGM")||names[k].equals("Cerebrum-GM")||names[k].equals("CerebralGM")||names[k].equals("Cerebellum-GM")||names[k].equals("GM")||names[k].equals("CerebellarGM"))
            for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)if(mask[x][y][z]){out[x][y][z][0]+=mem[x][y][z][k];if(seg[x][y][z]==4)out[x][y][z][2]=1.0f;}
        for(int z=0;z<Z;z++){
            for(int x=0;x<X;x++){out[x][0][z][2]=0;out[x][Y-1][z][2]=0;}
            for(int y=0;y<Y;y++){out[0][y][z][2]=0;out[X-1][y][z][2]=0;}
            boolean[][] seen=new boolean[X][Y];Queue<int[]> q=new ArrayDeque<>();q.add(new int[]{0,0});seen[0][0]=true;
            int[][] d={{1,0},{-1,0},{0,1},{0,-1}};
            while(!q.isEmpty()){int[] p=q.remove();for(int[] v:d){int xx=p[0]+v[0],yy=p[1]+v[1];if(xx>=0&&xx<X&&yy>=0&&yy<Y&&!seen[xx][yy]&&out[xx][yy][z][2]==0){seen[xx][yy]=true;q.add(new int[]{xx,yy});}}}
            for(int x=0;x<X;x++)for(int y=0;y<Y;y++)out[x][y][z][2]=seen[x][y]?0:1;
        }
        for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)for(int c=0;c<2;c++)out[x][y][z][c]*=255.0f;
        for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)out[x][y][z][2]*=255.0f;
    }
    static String bits(){StringBuilder s=new StringBuilder("[");boolean first=true;for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)for(int c=0;c<3;c++){if(!first)s.append(',');first=false;s.append(Integer.toUnsignedString(Float.floatToRawIntBits(out[x][y][z][c])));}return s.append(']').toString();}
    public static void main(String[] args){init();create();System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"data\":"+bits()+"}");}
}
