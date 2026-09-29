import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Formula excerpt for the six-face neighbor loop in propagateThinning. */
public final class ToadsThinningNeighborScheduleReferenceHarness {
    private static final int N=5,C=3,CX=2,CY=2,CZ=2,K=0;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final class Scenario {
        String name;byte[][][] seg=new byte[N][N][N];boolean[][][][] available=new boolean[N][N][N][C];
        float value,volume=101.0f,minimum=100.0f,limit=5.0f;List<int[]> seeds=new ArrayList<>();
        Scenario(String n){name=n;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)seg[x][y][z]=1;}
    }
    private static float mem(int x,int y,int z,int k){return (((x*5+y*7+z*11+k*3)%29)+1)/30.0f;}
    private static void allNeighbors(Scenario s,boolean setAvailable){for(int[] d:FACE){int x=CX+d[0],y=CY+d[1],z=CZ+d[2];s.seg[x][y][z]=K;s.available[x][y][z][K]=setAvailable;}}
    private static Scenario scenario(String name){
        Scenario s=new Scenario(name);
        if(name.equals("all_six")){allNeighbors(s,true);s.value=0.125f;}
        else if(name.equals("single_count_tie")){s.seg[1][2][2]=K;s.available[1][2][2][K]=true;s.value=-0.25f;s.limit=0.5f;s.seeds.add(new int[]{Float.floatToIntBits(-0.25f),4,4,4,2});}
        else if(name.equals("volume_equal")){allNeighbors(s,true);s.volume=s.minimum;s.value=0.125f;}
        else if(name.equals("limit_equal")){s.seg[1][2][2]=K;s.available[1][2][2][K]=true;s.value=-0.5f;s.limit=0.5f;}
        else if(name.equals("limit_above")){s.seg[1][2][2]=K;s.available[1][2][2][K]=true;s.value=-0.49f;s.limit=0.5f;}
        else if(name.equals("availability_gate")){allNeighbors(s,true);s.available[1][2][2][K]=false;s.value=0.125f;}
        return s;
    }
    private static String bits(float v){return Integer.toString(Float.floatToIntBits(v));}
    private static String labels(byte[][][] a){StringBuilder b=new StringBuilder("[");boolean first=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!first)b.append(',');first=false;b.append(a[x][y][z]);}return b.append(']').toString();}
    private static String available(boolean[][][][] a){StringBuilder b=new StringBuilder("[");boolean first=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){if(!first)b.append(',');first=false;b.append(a[x][y][z][k]);}return b.append(']').toString();}
    private static void run(Scenario s,boolean first){
        if(!first)System.out.print(',');
        BinaryHeap4D heap=new BinaryHeap4D(2,BinaryHeap4D.MAXTREE);
        for(int[] a:s.seeds)heap.addValue(Float.intBitsToFloat(a[0]),a[1],a[2],a[3],a[4]);
        boolean[][][][] before=new boolean[N][N][N][C];for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)before[x][y][z][k]=s.available[x][y][z][k];
        List<int[]> inserted=new ArrayList<>();
        for(int[] d:FACE){int xi=CX+d[0],yj=CY+d[1],zl=CZ+d[2];if(s.available[xi][yj][zl][K]){
            int count=0;for(int ni=-1;ni<=1;ni++)for(int nj=-1;nj<=1;nj++)for(int nl=-1;nl<=1;nl++)if(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]>0&&s.seg[xi+ni][yj+nj][zl+nl]==K)count++;
            float order;
            if(count==1)order=s.value;
            else {int b=0;float max=0.0f;for(int m=0;m<C;m++)if(mem(xi,yj,zl,m)>max){b=m;max=mem(xi,yj,zl,m);}float ratio=K==b?1.0f:(mem(xi,yj,zl,b)+1e-10f)/(mem(xi,yj,zl,K)+1e-10f);order=s.value-1.0f/ratio;}
            if(s.volume>s.minimum&&order>-s.limit){s.available[xi][yj][zl][K]=false;heap.addValue(order,xi,yj,zl,K);inserted.add(new int[]{Float.floatToIntBits(order),xi,yj,zl,K});}
        }}
        System.out.print("{\"name\":\""+s.name+"\",\"segmentation\":"+labels(new ScenarioCopy(s.seg).copy())+",\"available_before\":"+available(before)+",\"value_bits\":"+bits(s.value)+",\"volume_bits\":"+bits(s.volume)+",\"minimum_bits\":"+bits(s.minimum)+",\"limit_bits\":"+bits(s.limit)+",\"seeds\":[");
        for(int i=0;i<s.seeds.size();i++){if(i>0)System.out.print(',');int[] a=s.seeds.get(i);System.out.print("["+a[0]+","+a[1]+","+a[2]+","+a[3]+","+a[4]+"]");}
        System.out.print("],\"available_after\":"+available(s.available)+",\"inserted\":[");
        for(int i=0;i<inserted.size();i++){if(i>0)System.out.print(',');int[] a=inserted.get(i);System.out.print("["+a[0]+","+a[1]+","+a[2]+","+a[3]+","+a[4]+"]");}
        System.out.print("],\"popped\":[");boolean f=true;while(heap.isNotEmpty()){if(!f)System.out.print(',');f=false;System.out.print('['+bits(heap.getFirst())+','+heap.getFirstX()+','+heap.getFirstY()+','+heap.getFirstZ()+','+heap.getFirstK()+']');heap.removeFirst();}System.out.print("]}");
    }
    private static final class ScenarioCopy {byte[][][] a;ScenarioCopy(byte[][][] v){a=v;}byte[][][] copy(){byte[][][] b=new byte[N][N][N];for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)b[x][y][z]=a[x][y][z];return b;}}
    public static void main(String[] args){String[] names={"all_six","single_count_tie","volume_equal","limit_equal","limit_above","availability_gate"};System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.propagateThinning neighbor scheduling\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[5,5,5,3],\"cases\":[");for(int i=0;i<names.length;i++)run(scenario(names[i]),i==0);System.out.println("]}");}
}
