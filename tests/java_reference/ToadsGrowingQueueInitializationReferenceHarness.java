import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/**
 * Synthetic source excerpt for VandyCRUISE initGrowingLabels.
 * Selected MultipleToadSegmentation.java is dated April 2006 and names
 * Pierre-Louis Bazin as author.
 */
public final class ToadsGrowingQueueInitializationReferenceHarness {
    private static final int N=7,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] skeleton=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MINTREE);
    private final String kind;
    private ToadsGrowingQueueInitializationReferenceHarness(String kind){this.kind=kind;}
    private void initialize(){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            skeleton[x][y][z]=EMPTY;
            ordering[x][y][z]=-((x*13+y*7+z*3)%19)/10.0f;
            for(int k=0;k<C;k++)mems[x][y][z][k]=(((x*5+y*7+z*11+k*3)%29)+1)/30.0f;
        }
        if(kind.equals("mixed")){
            skeleton[3][3][2]=0;skeleton[3][3][4]=1;
            skeleton[2][3][3]=2;skeleton[4][3][3]=3;
        } else if(kind.equals("single_seed")) {
            skeleton[3][3][3]=0;
        } else if(kind.equals("no_boundary")) {
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)skeleton[x][y][z]=0;
        }
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int k=0;k<C;k++)available[x][y][z][k]=skeleton[x][y][z]==EMPTY;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){
            if(!available[x][y][z][k])continue;
            boolean boundary=false;
            for(int[] d:FACE)if(skeleton[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            available[x][y][z][k]=false;
            int best=0;float max=0.0f;
            for(int c=0;c<C;c++)if(mems[x][y][z][c]>max){best=c;max=mems[x][y][z][c];}
            float ratio=k==best?1.0f:(mems[x][y][z][best]+ZERO)/(mems[x][y][z][k]+ZERO);
            tree.addValue(ordering[x][y][z]+ratio,x,y,z,k);
        }
    }
    private static String labels(byte[][][] a){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(v==EMPTY?'a':(char)('0'+v));return b.toString();}
    private static String available(boolean[][][][] a){StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(a[x][y][z][k]?'1':'0');return b.toString();}
    private String queue(){List<String> out=new ArrayList<>();while(tree.isNotEmpty()){out.add("["+Float.floatToIntBits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}return "["+String.join(",",out)+"]";}
    private void emit(boolean first){if(!first)System.out.print(',');System.out.print("{\"name\":\""+kind+"\",\"skeleton\":\""+labels(skeleton)+"\",\"available\":\""+available(available)+"\",\"queue\":"+queue()+"}");}
    public static void main(String[] args){
        String[] cases={"mixed","single_seed","no_boundary"};
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.initGrowingLabels\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[7,7,7,4],\"cases\":[");
        for(int i=0;i<cases.length;i++){ToadsGrowingQueueInitializationReferenceHarness h=new ToadsGrowingQueueInitializationReferenceHarness(cases[i]);h.initialize();h.emit(i==0);}
        System.out.println("]}");
    }
}
