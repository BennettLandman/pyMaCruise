import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** One-event source excerpt for VandyCRUISE MultipleToadSegmentation.propagateGrowing. */
public final class ToadsGrowingQueueStepReferenceHarness {
    private static final int N=5,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final String REJECT="020033121011213000023100323";
    private final byte[][][] seg=new byte[N][N][N],skel=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final float[] volume={10,20,30,40};
    private final BinaryHeap4D tree=new BinaryHeap4D(8,BinaryHeap4D.MINTREE);
    private final String name;
    private int queuedLabel;
    private final List<String> inserted=new ArrayList<>();
    private ToadsGrowingQueueStepReferenceHarness(String name){this.name=name;}
    private void initialize(){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            seg[x][y][z]=0;skel[x][y][z]=0;
            for(int k=0;k<C;k++)mems[x][y][z][k]=(((x*5+y*7+z*11+k*3)%29)+1)/30.0f;
        }
        if(name.equals("direct_allowed")){
            seg[1][2][2]=skel[1][2][2]=1;skel[2][2][2]=EMPTY;queuedLabel=1;
            for(int[] d:FACE){int x=2+d[0],y=2+d[1],z=2+d[2];if(skel[x][y][z]!=1){skel[x][y][z]=EMPTY;available[x][y][z][1]=true;}}
            available[1][2][2][1]=false;
        }else if(name.equals("direct_rejected")){
            decode(REJECT,1,1,1);skel[2][2][2]=EMPTY;queuedLabel=3;
            for(int[] d:FACE){int x=2+d[0],y=2+d[1],z=2+d[2];if(seg[x][y][z]==3)available[x][y][z][3]=false;}
        }else {
            seg[2][2][2]=1;skel[2][2][2]=0;queuedLabel=3;
        }
        tree.addValue(.25f,2,2,2,queuedLabel);
        tree.addValue(.5f,3,3,3,2);
    }
    private void decode(String s,int ox,int oy,int oz){for(int i=0;i<27;i++){int x=i/9,y=(i/3)%3,z=i%3;seg[ox+x][oy+y][oz+z]=(byte)(s.charAt(i)-'0');skel[ox+x][oy+y][oz+z]=seg[ox+x][oy+y][oz+z];}}
    private short[][][] toShort(){short[][][] a=new short[N][N][N];for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)a[x][y][z]=skel[x][y][z];return a;}
    private boolean homeomorphic(int x,int y,int z,int k){short[][][] a=new short[N][N][N];for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)a[i][j][l]=seg[i][j][l];return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(a,(short)x,(short)y,(short)z,(short)k);}
    private Boolean topologyAccepted;
    private float poppedValue;
    private int poppedX,poppedY,poppedZ,poppedLabel;
    private void process(){
        float val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
        poppedValue=val;poppedX=x;poppedY=y;poppedZ=z;poppedLabel=k;
        if(skel[x][y][z]==EMPTY){
            boolean regular=homeomorphic(x,y,z,k);
            regular=homeomorphic(x,y,z,k); // Duplicate LUT call retained from selected source.
            topologyAccepted=regular;
            if(regular){
                skel[x][y][z]=(byte)k;for(int m=0;m<C;m++)available[x][y][z][m]=false;
                ordering[x][y][z]=val;volume[k]+=mems[x][y][z][k];seg[x][y][z]=(byte)k;
                for(int[] d:FACE){int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;
                    int best=0;float max=0;for(int c=0;c<C;c++)if(mems[xi][yj][zl][c]>max){best=c;max=mems[xi][yj][zl][c];}
                    float ratio=k==best?1.0f:(mems[xi][yj][zl][best]+ZERO)/(mems[xi][yj][zl][k]+ZERO);
                    float order=val+ratio;available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);
                    inserted.add("["+Float.floatToIntBits(order)+","+xi+","+yj+","+zl+","+k+"]");
                }
            }else available[x][y][z][k]=true;
        }
    }
    private static String labels(byte[][][] a){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append((char)('0'+v));return b.toString();}
    private static String skeleton(byte[][][] a){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(v==EMPTY?'a':(char)('0'+v));return b.toString();}
    private static String available(boolean[][][][] a){StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(a[x][y][z][k]?'1':'0');return b.toString();}
    private static String floats(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();}
    private String queue(){List<String> items=new ArrayList<>();while(tree.isNotEmpty()){items.add("["+Float.floatToIntBits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}return "["+String.join(",",items)+"]";}
    private void emit(boolean first){process();if(!first)System.out.print(',');System.out.print("{\"name\":\""+name+"\",\"topology_accepted\":"+(topologyAccepted==null?"null":topologyAccepted.toString())+",\"popped\":["+Float.floatToIntBits(poppedValue)+","+poppedX+","+poppedY+","+poppedZ+","+poppedLabel+"],\"segmentation\":\""+labels(seg)+"\",\"skeleton\":\""+skeleton(skel)+"\",\"available\":\""+available(available)+"\",\"ordering_bits\":"+floats(centerOrdering())+",\"volumes_bits\":"+floats(volume)+",\"inserted\":["+String.join(",",inserted)+"],\"queue_after\":"+queue()+"}");}
    private float[] centerOrdering(){return new float[]{ordering[2][2][2]};}
    public static void main(String[] args){ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();String[] names={"direct_allowed","direct_rejected","stale_skeleton"};System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.propagateGrowing one queue item\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"cases\":[");for(int i=0;i<names.length;i++){ToadsGrowingQueueStepReferenceHarness h=new ToadsGrowingQueueStepReferenceHarness(names[i]);h.initialize();h.emit(i==0);}System.out.println("]}");}
}
