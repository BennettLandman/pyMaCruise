import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected initThinningLabels + minVolume + propagateThinning synthetic pass. */
public final class LesionToadsThinningPassReferenceHarness {
    private static final int N=11,C=3,LESION=1;
    private static final float RATIO=0.4f,LIMIT=1000.0f,ZERO=1.0e-10f,EMPTY=-1.0f;
    private static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static short[][][] seg=new short[N][N][N],skeleton=new short[N][N][N];
    private static float[][][][] mem=new float[N][N][N][C],shape=new float[N][N][N][C];
    private static float[][][] lesion=new float[N][N][N],black=new float[N][N][N],ordering=new float[N][N][N];
    private static boolean[][][][] available=new boolean[N][N][N][C];
    private static float[] volume=new float[C],minimum=new float[C];

    private static int classAt(int x){return x<=4?0:(x<=7?1:2);}
    private static float score(int x,int y,int z,int k){return k==LESION?((mem[x][y][z][k]+lesion[x][y][z])+black[x][y][z]):mem[x][y][z][k];}
    private static int best(int x,int y,int z){int id=0;float max=0.0f;for(int k=0;k<C;k++){float v=score(x,y,z,k);if(v>max){id=k;max=mem[x][y][z][k];}}return id;}
    private static float ratio(int x,int y,int z,int current,int best){if(current==best)return 1.0f;return (score(x,y,z,best)+ZERO)/(score(x,y,z,current)+ZERO);}
    private static boolean boundary(int x,int y,int z,int k){return seg[x-1][y][z]!=k||seg[x+1][y][z]!=k||seg[x][y-1][z]!=k||seg[x][y+1][z]!=k||seg[x][y][z-1]!=k||seg[x][y][z+1]!=k;}
    private static String bits(float v){return Integer.toString(Float.floatToRawIntBits(v));}
    private static String shortArray(short[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(a[x][y][z]);}return b.append(']').toString();}
    private static String boolArray(){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){if(!f)b.append(',');f=false;b.append(available[x][y][z][k]);}return b.append(']').toString();}
    private static String floatBitsArray(float[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(Float.floatToRawIntBits(a[x][y][z]));}return b.append(']').toString();}
    public static void main(String[] args){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int label=classAt(x);
            seg[x][y][z]=(short)label;
            skeleton[x][y][z]=(short)(x==0||x==N-1||y==0||y==N-1||z==0||z==N-1?7:-1);
            lesion[x][y][z]=0.1f;black[x][y][z]=0.05f;ordering[x][y][z]=0.375f;
            for(int k=0;k<C;k++){
                mem[x][y][z][k]=label==k?0.8f:(k==1?0.1f:0.05f);
                shape[x][y][z][k]=0.0f;
                available[x][y][z][k]=((x+y+z+k)%2)==0;
            }
        }
        // Selected initThinningLabels: reset volumes, then scan the interior.
        for(int k=0;k<C;k++)volume[k]=0.0f;
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            int current=seg[x][y][z];skeleton[x][y][z]=(short)current;
            for(int k=0;k<C;k++){
                if(current==k){float inc=k==LESION?(mem[x][y][z][k]+lesion[x][y][z])+black[x][y][z]:mem[x][y][z][k];volume[k]+=inc;}
                available[x][y][z][k]=current==k;
            }
        }
        float[] initialVolume=volume.clone();
        BinaryHeap4D tree=new BinaryHeap4D(N*N*N*C,BinaryHeap4D.MAXTREE);
        List<String> inserted=new ArrayList<>();
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            int k=seg[x][y][z];if(boundary(x,y,z,k)){
                available[x][y][z][k]=false;float p=-1.0f/ratio(x,y,z,k,best(x,y,z));tree.addValue(p,x,y,z,k);
                inserted.add("["+bits(p)+","+x+","+y+","+z+","+k+"]");
            }
        }
        for(int k=0;k<C;k++)minimum[k]=Numerics.max(100.0f,RATIO*volume[k]);
        List<String> events=new ArrayList<>();float val=0.0f;
        while(tree.isNotEmpty()&&val>-LIMIT){
            short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();val=tree.getFirst();tree.removeFirst();boolean processed=false;
            if(volume[k]>minimum[k]&&skeleton[x][y][z]==k){
                int id=best(x,y,z);boolean changed=false;
                if(id==k)changed=false;
                else throw new AssertionError("synthetic pass should retain the current best label");
                available[x][y][z][k]=false;ordering[x][y][z]=val;
                volume[k]-=k==LESION?((mem[x][y][z][k]+lesion[x][y][z])+black[x][y][z]):mem[x][y][z][k];
                skeleton[x][y][z]=(short)EMPTY;if(changed)seg[x][y][z]=(short)id;processed=true;
                for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(available[xi][yj][zl][k]){float order=val-1.0f/ratio(xi,yj,zl,k,best(xi,yj,zl));if(volume[k]>minimum[k]&&order>-LIMIT){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);}}}
            }
            events.add("["+bits(val)+","+x+","+y+","+z+","+k+","+processed+"]");
        }
        StringBuilder ins=new StringBuilder("[");for(int i=0;i<inserted.size();i++){if(i>0)ins.append(',');ins.append(inserted.get(i));}ins.append(']');
        StringBuilder trace=new StringBuilder("[");for(int i=0;i<events.size();i++){if(i>0)trace.append(',');trace.append(events.get(i));}trace.append(']');
        StringBuilder rem=new StringBuilder("[");boolean first=true;while(tree.isNotEmpty()){if(!first)rem.append(',');first=false;rem.append("["+bits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}rem.append(']');
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"initial_volume_bits\":["+bits(initialVolume[0])+","+bits(initialVolume[1])+","+bits(initialVolume[2]) +"],\"initial_inserted\":"+ins+",\"minimum_bits\":["+bits(minimum[0])+","+bits(minimum[1])+","+bits(minimum[2]) +"],\"events\":"+trace+",\"remaining\":"+rem+",\"segmentation\":"+shortArray(seg)+",\"skeleton\":"+shortArray(skeleton)+",\"available\":"+boolArray()+",\"ordering_bits\":"+floatBitsArray(ordering)+",\"volume_bits\":["+bits(volume[0])+","+bits(volume[1])+","+bits(volume[2])+"]}");
    }
}
