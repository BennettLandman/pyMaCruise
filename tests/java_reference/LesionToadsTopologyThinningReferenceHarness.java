import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected LesionTOADS thinning excerpt exercised with topology-changing memberships. */
public final class LesionToadsTopologyThinningReferenceHarness {
    private static final int N=11, C=3, EMPTY=-1;
    private static final float RATIO=0.4f, LIMIT=10.0f, ZERO=1.0e-10f;
    private static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final short[][][] seg=new short[N][N][N], skel=new short[N][N][N];
    private static final float[][][][] mem=new float[N][N][N][C], shape=new float[N][N][N][C];
    private static final float[][][] lesion=new float[N][N][N], black=new float[N][N][N], order=new float[N][N][N];
    private static final boolean[][][][] available=new boolean[N][N][N][C];
    private static final float[] volume=new float[C], minimum=new float[C];

    private static int labelAt(int x,int y,int z,int state) {
        int base=((x/3)+2*(y/3)+(z/3))%3;
        state=(1103515245*state+12345)&0x7fffffff;
        return state%100<5?(state/100)%3:base;
    }
    private static boolean boundary(int x,int y,int z,int k) {
        return seg[x-1][y][z]!=k||seg[x+1][y][z]!=k||seg[x][y-1][z]!=k
            ||seg[x][y+1][z]!=k||seg[x][y][z-1]!=k||seg[x][y][z+1]!=k;
    }
    private static float ratio(int x,int y,int z,int current,int best) {
        if(current==best)return 1.0f;
        return (mem[x][y][z][best]+ZERO)/(mem[x][y][z][current]+ZERO);
    }
    private static int best(int x,int y,int z) {
        int id=0;float max=0.0f;
        for(int k=0;k<C;k++)if(mem[x][y][z][k]>max){id=k;max=mem[x][y][z][k];}
        return id;
    }
    private static String bits(float v){return Integer.toString(Float.floatToRawIntBits(v));}
    private static String shorts(short[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(a[x][y][z]);}return b.append(']').toString();}
    private static String bools(){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){if(!f)b.append(',');f=false;b.append(available[x][y][z][k]);}return b.append(']').toString();}
    private static String floatBits(float[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){if(!f)b.append(',');f=false;b.append(Float.floatToRawIntBits(a[x][y][z]));}return b.append(']').toString();}

    public static void main(String[] args) {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        int state=1;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int k=labelAt(x,y,z,state);state=(1103515245*state+12345)&0x7fffffff;
            seg[x][y][z]=(short)k;
            boolean border=x==0||x==N-1||y==0||y==N-1||z==0||z==N-1;
            skel[x][y][z]=(short)(border?7:EMPTY);
            for(int m=0;m<C;m++)mem[x][y][z][m]=0.005f;
            mem[x][y][z][k]=0.49f;mem[x][y][z][(k+1)%C]=0.51f;
        }
        for(int k=0;k<C;k++)volume[k]=0.0f;
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            int current=seg[x][y][z];skel[x][y][z]=(short)current;
            for(int k=0;k<C;k++){if(current==k)volume[k]+=mem[x][y][z][k];available[x][y][z][k]=current==k;}
        }
        float[] initialVolume=volume.clone();
        BinaryHeap4D tree=new BinaryHeap4D(N*N*N*C,BinaryHeap4D.MAXTREE);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){
            int k=seg[x][y][z];if(boundary(x,y,z,k)){available[x][y][z][k]=false;tree.addValue(-1.0f/ratio(x,y,z,k,best(x,y,z)),x,y,z,k);}
        }
        for(int k=0;k<C;k++)minimum[k]=Numerics.max(100.0f,RATIO*volume[k]);
        List<String> events=new ArrayList<>();int accepted=0,rejected=0,changes=0,insertions=0;float val=0.0f;
        while(tree.isNotEmpty()&&val>-LIMIT){
            short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();val=tree.getFirst();tree.removeFirst();
            boolean processed=false;int candidate=-1;Boolean topology=null;List<String> added=new ArrayList<>();
            if(volume[k]>minimum[k]&&skel[x][y][z]==k){
                int id=best(x,y,z);boolean direct=false;int next=-1;float prior=0.0f;
                if(id!=k){
                    for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l==1){
                        short neighbor=seg[x+i][y+j][z+l];if(neighbor==id)direct=true;
                        else if(neighbor!=k&&shape[x+i][y+j][z+l][id]>prior){next=neighbor;prior=shape[x+i][y+j][z+l][id];}
                    }
                    if(direct)candidate=id;else if(next!=-1)candidate=next;
                    if(candidate!=-1)topology=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(seg,x,y,z,(short)candidate);
                    if(Boolean.TRUE.equals(topology)){seg[x][y][z]=(short)candidate;accepted++;changes++;}
                    else if(Boolean.FALSE.equals(topology))rejected++;
                }
                available[x][y][z][k]=false;order[x][y][z]=val;volume[k]-=mem[x][y][z][k];skel[x][y][z]=(short)EMPTY;processed=true;
                for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){
                    int xi=x+i,yj=y+j,zl=z+l;
                    if(available[xi][yj][zl][k]){float priority=val-1.0f/ratio(xi,yj,zl,k,best(xi,yj,zl));
                        if(volume[k]>minimum[k]&&priority>-LIMIT){available[xi][yj][zl][k]=false;tree.addValue(priority,xi,yj,zl,k);added.add("["+bits(priority)+","+xi+","+yj+","+zl+","+k+"]");insertions++;}
                    }
                }
            }
            String top=topology==null?"null":topology.toString();String cand=candidate<0?"null":Integer.toString(candidate);
            events.add("["+bits(val)+","+x+","+y+","+z+","+k+","+processed+","+cand+","+top+","+added+"]");
        }
        StringBuilder rem=new StringBuilder("[");boolean first=true;while(tree.isNotEmpty()){if(!first)rem.append(',');first=false;rem.append("["+bits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}rem.append(']');
        StringBuilder ev=new StringBuilder("[");for(int i=0;i<events.size();i++){if(i>0)ev.append(',');ev.append(events.get(i));}ev.append(']');
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"initial_volume_bits\":["+bits(initialVolume[0])+","+bits(initialVolume[1])+","+bits(initialVolume[2])+
            "],\"minimum_bits\":["+bits(minimum[0])+","+bits(minimum[1])+","+bits(minimum[2])+
            "],\"events\":"+ev+",\"remaining\":"+rem+",\"segmentation\":"+shorts(seg)+",\"skeleton\":"+shorts(skel)+
            ",\"available\":"+bools()+",\"ordering_bits\":"+floatBits(order)+",\"volume_bits\":["+bits(volume[0])+","+bits(volume[1])+","+bits(volume[2])+
            "],\"topology_accepted_count\":"+accepted+",\"topology_rejected_count\":"+rejected+",\"changed_label_count\":"+changes+",\"inserted_neighbor_count\":"+insertions+"}");
    }
}
