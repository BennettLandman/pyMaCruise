import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/**
 * Source-order composition reference: run the existing selected growth
 * harness, refresh relation adjacency, then initialize/run the next thinning
 * queue with the short-wrapper defaults. Registration and shape transforms
 * are fixed; zero shape priors are supplied at this seam.
 */
public final class ToadsDefaultContinueThinningReferenceHarness {
    private static final int N=9,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f,LIMIT=3.0f,RATIO=.5f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] seg=new byte[N][N][N],skel=new byte[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final boolean[][] relations=new boolean[C][C];
    private final boolean[][][] relationAt=new boolean[N][N][N];
    private final float[][][][] shape=new float[N][N][N][C];
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MAXTREE);
    private final float[] volume=new float[C],minimum=new float[C],initial=new float[C];
    private final int[] boundaries=new int[C];
    private final List<String> processed=new ArrayList<>();

    private static Object field(Object object,String name)throws Exception {
        Field f=object.getClass().getDeclaredField(name);f.setAccessible(true);return f.get(object);
    }
    private void initializeFromSelectedGrowth() throws Exception {
        ToadsGrowingPassReferenceHarness growth=new ToadsGrowingPassReferenceHarness();
        Method initialize=growth.getClass().getDeclaredMethod("initialize");initialize.setAccessible(true);initialize.invoke(growth);
        Method propagate=growth.getClass().getDeclaredMethod("propagate");propagate.setAccessible(true);propagate.invoke(growth);
        byte[][][] sourceSeg=(byte[][][])field(growth,"seg");float[][][][] sourceMems=(float[][][][])field(growth,"mems");
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            seg[x][y][z]=sourceSeg[x][y][z];
            for(int k=0;k<C;k++)mems[x][y][z][k]=sourceMems[x][y][z][k];
        }
        // Relation rows 0..C-1 are the selected buildRelationList class adjacency.
        boolean[][] fullRelation=ToadsRelationListReferenceHarness.build(seg,C).relations;
        for(int k=0;k<C;k++)System.arraycopy(fullRelation[k],0,relations[k],0,C);
        // The event kernel reads the refreshed relation matrix; the spatial
        // transformed-shape values are fixed zero arrays for this fixture.
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)skel[x][y][z]=(byte)EMPTY;
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int k=seg[x][y][z];skel[x][y][z]=(byte)k;volume[k]+=mems[x][y][z][k];
            for(int c=0;c<C;c++)available[x][y][z][c]=c==k;
        }
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int k=seg[x][y][z];boolean boundary=false;int neighborCount=0,ni=0,nj=0,nl=0;
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++) {
                int d2=i*i+j*j+l*l;
                if(d2==1&&seg[x+i][y+j][z+l]!=k)boundary=true;
                else if(d2>0){neighborCount++;ni=x+i;nj=y+j;nl=z+l;}
            }
            if(!boundary)continue;
            available[x][y][z][k]=false;
            if(neighborCount==1) {
                neighborCount=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(i*i+j*j+l*l>0&&seg[ni+i][nj+j][nl+l]==k)neighborCount++;
                if(neighborCount<=2)throw new AssertionError("selected source endpoint precondition became reachable");
            }
            int best=0;float max=0;for(int c=0;c<C;c++)if(mems[x][y][z][c]>max){best=c;max=mems[x][y][z][c];}
            float ratio=k==best?1.0f:(mems[x][y][z][best]+ZERO)/(mems[x][y][z][k]+ZERO);
            tree.addValue(-1.0f/ratio,x,y,z,k);boundaries[k]++;
        }
        for(int k=0;k<C;k++){minimum[k]=Math.max(100.0f,RATIO*volume[k]);initial[k]=volume[k];}
    }
    private int bestMembership(int x,int y,int z){int best=0;float max=0;for(int k=0;k<C;k++)if(mems[x][y][z][k]>max){best=k;max=mems[x][y][z][k];}return best;}
    private boolean topology(int x,int y,int z,int target){
        short[][][] copy=new short[N][N][N];for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)copy[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(copy,(short)x,(short)y,(short)z,(short)target);
    }
    private void oneEvent(){
        float val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();
        processed.add("["+Float.floatToIntBits(val)+","+x+","+y+","+z+","+k+"]");tree.removeFirst();
        if(volume[k]>minimum[k]&&skel[x][y][z]==k){
            int id=bestMembership(x,y,z),next=EMPTY;boolean direct=false;float bestPrior=-.1f;
            if(id!=k)for(int[] d:FACE){
                int xi=x+d[0],yj=y+d[1],zl=z+d[2],neighbor=seg[xi][yj][zl];
                if(neighbor==id)direct=true;
                else if(neighbor!=k&&relations[neighbor][id]&&shape[xi][yj][zl][id]>bestPrior){next=neighbor;bestPrior=shape[xi][yj][zl][id];}
            }
            boolean change=id!=k&&(direct||next!=EMPTY);
            if(change){if(!direct)id=next;change=topology(x,y,z,id);}
            available[x][y][z][k]=false;volume[k]-=mems[x][y][z][k];skel[x][y][z]=(byte)EMPTY;
            if(change)seg[x][y][z]=(byte)id;
            for(int[] d:FACE){
                int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;
                int count=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]>0&&seg[xi+i][yj+j][zl+l]==k)count++;
                float order;
                if(count==1)order=val;
                else {int best=bestMembership(xi,yj,zl);float ratio=k==best?1.0f:(mems[xi][yj][zl][best]+ZERO)/(mems[xi][yj][zl][k]+ZERO);order=val-1.0f/ratio;}
                if(volume[k]>minimum[k]&&order>-LIMIT){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);}
            }
        }
    }
    private void run(){
        float val=0.0f;
        while(tree.isNotEmpty()&&val>-LIMIT){
            val=tree.getFirst();oneEvent();
        }
    }
    private String queue(){List<String> q=new ArrayList<>();while(tree.isNotEmpty()){
        q.add("["+Float.floatToIntBits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}
        return "["+String.join(",",q)+"]";}
    private static String labels(byte[][][] a,boolean empty){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(empty&&v==EMPTY?'a':(char)('0'+v));return b.toString();}
    private static String available(boolean[][][][] a){StringBuilder b=new StringBuilder();for(boolean[][][] x:a)for(boolean[][] y:x)for(boolean[] z:y)for(boolean v:z)b.append(v?'1':'0');return b.toString();}
    private static String bools(boolean[][] a){StringBuilder b=new StringBuilder("[");boolean first=true;for(boolean[] x:a)for(boolean v:x){if(!first)b.append(',');first=false;b.append(v?'1':'0');}return b.append(']').toString();}
    private static String floats(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();}
    public static void main(String[] args)throws Exception {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        ToadsDefaultContinueThinningReferenceHarness h=new ToadsDefaultContinueThinningReferenceHarness();h.initializeFromSelectedGrowth();h.run();
        System.out.print("{\"source_order\":[\"continue.growing.segmentation\",\"continue.relation_state.class_adjacency\",\"exact_memberships\",\"fixed_transformed_shapes\",\"initThinningLabels\",\"propagateThinning\"],\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"shape\":["+N+","+N+","+N+","+C+"],\"volume_ratio_bits\":"+Float.floatToIntBits(RATIO)+",\"limit_bits\":"+Float.floatToIntBits(LIMIT)+",\"relation_adjacency_bits\":\""+bools(h.relations)+"\",\"initial_volumes_bits\":"+floats(h.initial)+",\"minimum_volumes_bits\":"+floats(h.minimum)+",\"boundary_counts\":["+h.boundaries[0]+","+h.boundaries[1]+","+h.boundaries[2]+","+h.boundaries[3]+"],\"processed_count\":"+h.processed.size()+",\"processed\":["+String.join(",",h.processed)+"],\"segmentation\":\""+labels(h.seg,false)+"\",\"skeleton\":\""+labels(h.skel,true)+"\",\"available\":\""+available(h.available)+"\",\"final_volumes_bits\":"+floats(h.volume)+",\"queue_after\":"+h.queue()+"}");
    }
}
