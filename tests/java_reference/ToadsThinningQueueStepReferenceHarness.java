import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/**
 * Executes one source-ordered propagateThinning event on synthetic grids.
 * Formula order is taken from VandyCRUISE MultipleToadSegmentation.java,
 * version April 2006, author Pierre-Louis Bazin, method propagateThinning.
 * This test harness uses the selected source CriticalPointLUT and BinaryHeap4D.
 */
public final class ToadsThinningQueueStepReferenceHarness {
    private static final int N=5, C=4, EMPTY=-1;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final String REJECT="020033121011213000023100323";

    private static final class State {
        int[][][] seg=new int[N][N][N], skel=new int[N][N][N];
        float[][][][] mem=new float[N][N][N][C], shape=new float[N][N][N][C];
        boolean[][][][] available=new boolean[N][N][N][C];
        boolean[][] relation=new boolean[C][C];
        float[] volume={200,210,220,230}, min={100,100,100,100};
        BinaryHeap4D queue=new BinaryHeap4D(8,BinaryHeap4D.MAXTREE);
        String name;
        State(String name) {
            this.name=name;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)
                mem[x][y][z]=new float[]{.1f,.8f,.05f,.05f};
            queue.addValue(-.25f,2,2,2,0);
            queue.addValue(-.5f,3,3,3,2);
            if(name.equals("direct_allowed")) {
                seg[1][2][2]=skel[1][2][2]=1;
            } else if(name.equals("direct_rejected")) {
                decode(REJECT,1,1,1);
                mem[2][2][2]=new float[]{.1f,.1f,.1f,.8f};
                queue.reset(); queue.addValue(-.25f,2,2,2,1); queue.addValue(-.5f,3,3,3,2);
            } else if(name.equals("related_allowed")) {
                seg[1][2][2]=skel[1][2][2]=2;
                relation[2][1]=true; shape[1][2][2][1]=.7f;
            } else if(name.equals("no_connection")) {
                seg[1][2][2]=skel[1][2][2]=2;
            } else if(name.equals("below_floor")) {
                volume[0]=100;
                seg[1][2][2]=skel[1][2][2]=1;
            } else if(name.equals("stale_skeleton")) {
                skel[2][2][2]=1;
                seg[1][2][2]=skel[1][2][2]=1;
            } else if(name.equals("queued_label_differs_from_segmentation")) {
                seg[2][2][2]=2; skel[2][2][2]=0;
            }
            int k=name.equals("direct_rejected")?1:0;
            for(int[] d:FACE)available[2+d[0]][2+d[1]][2+d[2]][k]=true;
        }
        void decode(String text,int ox,int oy,int oz) {
            for(int i=0;i<27;i++) {
                int x=i/9,y=(i/3)%3,z=i%3;
                seg[ox+x][oy+y][oz+z]=text.charAt(i)-'0';
                skel[ox+x][oy+y][oz+z]=seg[ox+x][oy+y][oz+z];
            }
        }
    }

    private static void process(State s) {
        float val=s.queue.getFirst();
        int x=s.queue.getFirstX(),y=s.queue.getFirstY(),z=s.queue.getFirstZ(),k=s.queue.getFirstK();
        s.queue.removeFirst();
        if(s.volume[k]>s.min[k] && s.skel[x][y][z]==k) {
            int id=0; float max=0;
            for(int c=0;c<C;c++)if(s.mem[x][y][z][c]>max){id=c;max=s.mem[x][y][z][c];}
            boolean direct=false,change=false; int next=EMPTY; float bestPrior=-.1f;
            if(id!=k)for(int[] d:FACE) {
                int xi=x+d[0],yj=y+d[1],zl=z+d[2],neighbor=s.seg[xi][yj][zl];
                if(neighbor==id)direct=true;
                else if(neighbor!=k && s.relation[neighbor][id] && s.shape[xi][yj][zl][id]>bestPrior) {
                    next=neighbor; bestPrior=s.shape[xi][yj][zl][id];
                }
            }
            if(id!=k && (direct || next!=EMPTY)) {
                change=true; if(!direct)id=next;
                change=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(
                    toShort(s.seg),(short)x,(short)y,(short)z,(short)id);
            }
            s.available[x][y][z][k]=false;
            s.volume[k]-=s.mem[x][y][z][k];
            s.skel[x][y][z]=EMPTY;
            if(change)s.seg[x][y][z]=id;
            for(int[] d:FACE) {
                int xi=x+d[0],yj=y+d[1],zl=z+d[2];
                if(!s.available[xi][yj][zl][k])continue;
                int count=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]>0 && s.seg[xi+i][yj+j][zl+l]==k)count++;
                float order;
                if(count==1)order=val;
                else {
                    int b=0;max=0;
                    for(int c=0;c<C;c++)if(s.mem[xi][yj][zl][c]>max){b=c;max=s.mem[xi][yj][zl][c];}
                    float ratio=k==b?1.0f:(s.mem[xi][yj][zl][b]+1e-10f)/(s.mem[xi][yj][zl][k]+1e-10f);
                    order=val-1.0f/ratio;
                }
                if(s.volume[k]>s.min[k] && order>-100f) {
                    s.available[xi][yj][zl][k]=false; s.queue.addValue(order,xi,yj,zl,k);
                }
            }
        }
    }
    private static short[][][] toShort(int[][][] a) {
        short[][][] b=new short[N][N][N];
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)b[x][y][z]=(short)a[x][y][z];
        return b;
    }
    private static String labels(int[][][] a) {
        StringBuilder b=new StringBuilder();for(int[][] x:a)for(int[] y:x)for(int v:y)b.append((char)('0'+v));return b.toString();
    }
    private static String skeleton(int[][][] a) {
        StringBuilder b=new StringBuilder();for(int[][] x:a)for(int[] y:x)for(int v:y)b.append(v<0?'a':(char)('0'+v));return b.toString();
    }
    private static String availability(boolean[][][][] a) {
        StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(a[x][y][z][k]?'1':'0');return b.toString();
    }
    private static String floats(float[] a) {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.toString(a[i]));}return b.append(']').toString();
    }
    private static String queue(BinaryHeap4D q) {
        List<String> items=new ArrayList<>();
        while(q.isNotEmpty()) {
            int bits=Float.floatToIntBits(q.getFirst());
            items.add("["+bits+","+q.getFirstX()+","+q.getFirstY()+","+q.getFirstZ()+","+q.getFirstK()+"]");q.removeFirst();
        }
        return "["+String.join(",",items)+"]";
    }
    private static void emit(State s,boolean first) {
        String inSeg=labels(s.seg),inSkel=skeleton(s.skel),inAvailable=availability(s.available);
        float[] inVolumes=s.volume.clone(); int k=s.name.equals("direct_rejected")?1:0;
        int priority=Float.floatToIntBits(-.25f);float[] center=s.mem[2][2][2].clone();
        process(s);
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+s.name+"\",\"input\":{\"segmentation\":\""+inSeg+"\",\"skeleton\":\""+inSkel+
            "\",\"available\":\""+inAvailable+"\",\"volumes\":"+floats(inVolumes)+",\"center_membership\":"+floats(center)+
            ",\"queue\":[["+priority+",2,2,2,"+k+"],["+Float.floatToIntBits(-.5f)+",3,3,3,2]],\"minimums\":[100,100,100,100],\"limit\":100},\"segmentation\":\""+
            labels(s.seg)+"\",\"skeleton\":\""+skeleton(s.skel)+"\",\"available\":\""+availability(s.available)+"\",\"volumes\":"+
            floats(s.volume)+",\"queue_after\":"+queue(s.queue)+"}");
    }
    public static void main(String[] args) {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        String[] names={"direct_allowed","direct_rejected","related_allowed","no_connection","below_floor","stale_skeleton","queued_label_differs_from_segmentation"};
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.propagateThinning one queue event\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"lut_source_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"cases\":[");
        for(int i=0;i<names.length;i++)emit(new State(names[i]),i==0);
        System.out.println("]}");
    }
}
