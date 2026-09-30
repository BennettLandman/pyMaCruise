import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Selected six-face post-removal loop from LesionToadSegmentation.propagateThinning. */
public final class LesionToadsThinningScheduleReferenceHarness {
    private static final int N=5,C=3,X=2,Y=2,Z=2;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";

    private static final class Scenario {
        final String name; final int removedClass; final boolean useMax;
        final float[][][][] mem=new float[N][N][N][C];
        final float[][][] lesion=new float[N][N][N], black=new float[N][N][N];
        final boolean[][][][] available=new boolean[N][N][N][C];
        float value=-0.25f, volume=8.0f, minimum=2.0f, limit=100.0f;
        Scenario(String name,int k,boolean max) {
            this.name=name; removedClass=k; useMax=max;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
                lesion[x][y][z]=(1+((x*3+y*5+z*7)%9))/32.0f;
                black[x][y][z]=(1+((x*5+y*3+z*2)%7))/64.0f;
                for(int c=0;c<C;c++)mem[x][y][z][c]=(2+((x*7+y*5+z*3+c*11)%15))/32.0f;
            }
            for(int[] d:FACE)available[X+d[0]][Y+d[1]][Z+d[2]][k]=true;
            if(name.equals("ratio_tie_limit_equal")||name.equals("ratio_tie_limit_above")){
                value=0.0f;limit=name.endsWith("equal")?1.0f:1.01f;
                for(int[] d:FACE){int x=X+d[0],y=Y+d[1],z=Z+d[2];for(int c=0;c<C;c++)mem[x][y][z][c]=0.25f;lesion[x][y][z]=0.0f;black[x][y][z]=0.0f;}
            }
        if(name.equals("lesion_max")||name.equals("lesion_sum")){value=-0.125f;volume=8.5f;minimum=2.25f;}
            if(name.equals("volume_equal")){volume=minimum=4.0f;}
            if(name.equals("volume_below")){volume=3.99f;minimum=4.0f;}
            if(name.equals("availability_gate"))available[X-1][Y][Z][k]=false;
        }
    }
    private static float score(Scenario s,int x,int y,int z,int c){
        if(c!=1&&c!=2)return s.mem[x][y][z][c];
        if(s.useMax){float b=s.lesion[x][y][z]>s.black[x][y][z]?s.lesion[x][y][z]:s.black[x][y][z];return s.mem[x][y][z][c]>b?s.mem[x][y][z][c]:b;}
        return (s.mem[x][y][z][c]+s.lesion[x][y][z])+s.black[x][y][z];
    }
    private static int best(Scenario s,int x,int y,int z){int best=0;float max=0.0f;for(int c=0;c<C;c++){if(score(s,x,y,z,c)>max){best=c;max=s.mem[x][y][z][c];}}return best;}
    private static float ratio(Scenario s,int x,int y,int z,int k,int b){if(k==b)return 1.0f;return (score(s,x,y,z,b)+1e-10f)/(score(s,x,y,z,k)+1e-10f);}
    private static String bits(float f){return Integer.toString(Float.floatToRawIntBits(f));}
    private static void run(Scenario s,boolean first){
        if(!first)System.out.print(',');
        BinaryHeap4D heap=new BinaryHeap4D(16,BinaryHeap4D.MAXTREE);
        List<int[]> inserted=new ArrayList<>();
        boolean[][][][] after=new boolean[N][N][N][C];
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int c=0;c<C;c++)after[x][y][z][c]=s.available[x][y][z][c];
        for(int[] d:FACE){int x=X+d[0],y=Y+d[1],z=Z+d[2],k=s.removedClass;if(!after[x][y][z][k])continue;
            float r=ratio(s,x,y,z,k,best(s,x,y,z));float order=s.value-1.0f/r;
            if(s.volume>s.minimum&&order>-s.limit){after[x][y][z][k]=false;heap.addValue(order,x,y,z,k);inserted.add(new int[]{Float.floatToRawIntBits(order),x,y,z,k});}
        }
        System.out.print("{\"name\":\""+s.name+"\",\"value_bits\":"+bits(s.value)+",\"volume_bits\":"+bits(s.volume)+",\"minimum_bits\":"+bits(s.minimum)+",\"limit_bits\":"+bits(s.limit)+",\"available_after\":[");
        boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int c=0;c<C;c++){if(!f)System.out.print(',');f=false;System.out.print(after[x][y][z][c]);}
        System.out.print("],\"inserted\":[");for(int i=0;i<inserted.size();i++){if(i>0)System.out.print(',');int[] a=inserted.get(i);System.out.print("["+a[0]+","+a[1]+","+a[2]+","+a[3]+","+a[4]+"]");}
        System.out.print("],\"popped\":[");f=true;while(heap.isNotEmpty()){if(!f)System.out.print(',');f=false;System.out.print("["+bits(heap.getFirst())+","+heap.getFirstX()+","+heap.getFirstY()+","+heap.getFirstZ()+","+heap.getFirstK()+"]");heap.removeFirst();}System.out.print("]}");
    }
    public static void main(String[] args){Scenario[] cases={
        new Scenario("all_faces_ordinary",0,false),new Scenario("lesion_max",1,true),
        new Scenario("lesion_sum",1,false),
        new Scenario("ratio_tie_limit_equal",0,false),new Scenario("ratio_tie_limit_above",0,false),
        new Scenario("volume_equal",0,false),new Scenario("volume_below",0,false),
        new Scenario("availability_gate",0,false)};
        System.out.print("{\"source\":\"LesionToadSegmentation.propagateThinning neighbor loop\",\"source_sha256\":\""+SOURCE_SHA+"\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[5,5,5,3],\"cases\":[");for(int i=0;i<cases.length;i++)run(cases[i],i==0);System.out.println("]}");}
}
