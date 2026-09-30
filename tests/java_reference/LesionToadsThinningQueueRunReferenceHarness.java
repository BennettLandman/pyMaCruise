import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Repeated selected LesionToadSegmentation.propagateThinning outer-loop excerpt. */
public final class LesionToadsThinningQueueRunReferenceHarness {
    private static final int N=5,C=3;
    private static final float LIMIT=2.0f;
    private static short[][][] seg=new short[N][N][N],skeleton=new short[N][N][N];
    private static float[][][][] mem=new float[N][N][N][C],shape=new float[N][N][N][C];
    private static float[][][] lesion=new float[N][N][N],black=new float[N][N][N],ordering=new float[N][N][N];
    private static boolean[][][][] available=new boolean[N][N][N][C];
    private static float[] volume={20.0f,1.0f,0.0f},minimum={0.0f,0.0f,0.0f};
    private static float score(int x,int y,int z,int c){return mem[x][y][z][c];}
    private static int best(int x,int y,int z){int id=0;float max=0.0f;for(int c=0;c<C;c++)if(score(x,y,z,c)>max){id=c;max=mem[x][y][z][c];}return id;}
    private static float ratio(int x,int y,int z,int k,int b){return k==b?1.0f:(score(x,y,z,b)+1e-10f)/(score(x,y,z,k)+1e-10f);}
    private static String bits(float v){return Integer.toString(Float.floatToRawIntBits(v));}
    private static String bools(){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int c=0;c<C;c++){if(!f)b.append(',');f=false;b.append(available[x][y][z][c]);}return b.append(']').toString();}
    private static String shorts(short[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(a[x][y][z]);}return b.append(']').toString();}
    private static String floats(float[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(Float.floatToRawIntBits(a[x][y][z]));}return b.append(']').toString();}
    public static void main(String[] args){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            for(int c=0;c<C;c++){mem[x][y][z][c]=c==0?0.5f:(c==1?0.25f:0.125f);shape[x][y][z][c]=0.0f;}
            if(x>0&&x<N-1&&y>0&&y<N-1&&z>0&&z<N-1&&!(x==2&&y==2&&z==2)&&!(x==1&&y==1&&z==1)&&!(x==1&&y==3&&z==1))available[x][y][z][0]=true;
        }
        skeleton[3][1][3]=2;
        BinaryHeap4D tree=new BinaryHeap4D(256,BinaryHeap4D.MAXTREE);tree.addValue(0.5f,2,2,2,0);tree.addValue(-1.0f,3,3,3,1);tree.addValue(-1.25f,3,1,3,2);tree.addValue(-2.0f,1,1,1,0);tree.addValue(-2.5f,1,3,1,0);
        List<String> processed=new ArrayList<>();float val=0.0f;
        while(tree.isNotEmpty()&&val>-LIMIT){
            short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();val=tree.getFirst();tree.removeFirst();
            boolean did=false;
            if(volume[k]>minimum[k]&&skeleton[x][y][z]==k){
                int id=best(x,y,z);boolean change=false;
                // Every source voxel has best label 0 here; this fixture's
                // processing label is 0, so no topology operation is needed.
                if(id!=k)throw new AssertionError("fixture should keep its current best label");
                available[x][y][z][k]=false;ordering[x][y][z]=val;volume[k]-=mem[x][y][z][k];skeleton[x][y][z]=-1;did=true;
                for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(available[xi][yj][zl][k]){float order=val-1.0f/ratio(xi,yj,zl,k,best(xi,yj,zl));if(volume[k]>minimum[k]&&order>-LIMIT){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);}}}
            }
            processed.add("["+bits(val)+","+x+","+y+","+z+","+k+","+did+"]");
        }
        StringBuilder trace=new StringBuilder("[");for(int i=0;i<processed.size();i++){if(i>0)trace.append(',');trace.append(processed.get(i));}trace.append(']');
        StringBuilder remaining=new StringBuilder("[");boolean first=true;while(tree.isNotEmpty()){if(!first)remaining.append(',');first=false;remaining.append("["+bits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}remaining.append(']');
        System.out.print("{\"source\":\"LesionToadSegmentation.propagateThinning outer loop\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"limit_bits\":"+bits(LIMIT)+",\"processed\":"+trace+",\"remaining\":"+remaining+",\"segmentation\":"+shorts(seg)+",\"skeleton\":"+shorts(skeleton)+",\"available\":"+bools()+",\"ordering_bits\":"+floats(ordering)+",\"volume_bits\":["+bits(volume[0])+","+bits(volume[1])+","+bits(volume[2])+"]}");
    }
}
