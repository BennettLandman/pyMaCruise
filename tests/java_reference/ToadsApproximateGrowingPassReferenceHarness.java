import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Sparse selected-source useApprox initGrowingLabels + propagateGrowing fixture. */
public final class ToadsApproximateGrowingPassReferenceHarness {
    private static final int N=9,C=3,A=2,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] seg=new byte[N][N][N],skel=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][A];
    private final byte[][][][] best=new byte[N][N][N][A];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final float[] volume={1,2,3};
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MINTREE);
    private final List<String> events=new ArrayList<>();

    private static int cls(int x,int y,int z){
        int fold=(x*3+y*2)%3-1,dx=x-4,dy=y-4,dz=z-4,radius=2+fold;
        if(dx*dx+dy*dy+dz*dz<radius*radius)return 1;
        int foldZ=4+(x+2*y)%3-1;
        return Math.abs(z-foldZ)<=1&&Math.abs(y-4)<=3&&Math.abs(x-4)<=3?2:0;
    }
    private void initialize(){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            int current=cls(x,y,z);seg[x][y][z]=(byte)current;skel[x][y][z]=EMPTY;
            ordering[x][y][z]=-((x*13+y*7+z*3)%19)/10.0f;
            for(int a=0;a<A;a++){
                best[x][y][z][a]=(byte)(a==0?current:(current+1)%C);
                int numerator=(x*17+y*11+z*7+a*13)%39+1;
                mems[x][y][z][a]=numerator/100.0f;
            }
            if(x>0&&x<N-1&&y>0&&y<N-1&&z>0&&z<N-1&&(x*2+y*3+z*5)%7==0)
                skel[x][y][z]=(byte)current;
        }
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int k=0;k<C;k++)available[x][y][z][k]=skel[x][y][z]==EMPTY;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){
            if(!available[x][y][z][k])continue;
            boolean boundary=false;for(int[] d:FACE)if(skel[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            available[x][y][z][k]=false;
            float ratio;if(k==best[x][y][z][0])ratio=1.0f;else{
                float cur=0.0f;for(int a=1;a<A;a++)if(best[x][y][z][a]==k)cur=mems[x][y][z][a];
                ratio=(mems[x][y][z][0]+ZERO)/(cur+ZERO);
            }
            tree.addValue(ordering[x][y][z]+ratio,x,y,z,k);
        }
    }
    private boolean homeomorphic(int x,int y,int z,int k){
        short[][][] a=new short[N][N][N];
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)a[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(a,(short)x,(short)y,(short)z,(short)k);
    }
    private static String item(float v,int x,int y,int z,int k){return "["+Float.floatToIntBits(v)+","+x+","+y+","+z+","+k+"]";}
    private void propagate(){
        while(tree.isNotEmpty()){
            float val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
            Boolean regular=null;List<String> inserted=new ArrayList<>();
            if(skel[x][y][z]==EMPTY){
                regular=homeomorphic(x,y,z,k);
                if(regular){
                    skel[x][y][z]=(byte)k;for(int m=0;m<C;m++)available[x][y][z][m]=false;
                    ordering[x][y][z]=val;
                    for(int a=0;a<A;a++)if(best[x][y][z][a]==k)volume[k]+=mems[x][y][z][a];
                    seg[x][y][z]=(byte)k;
                    for(int[] d:FACE){int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;
                        float ratio;if(k==best[xi][yj][zl][0])ratio=1.0f;else{
                            float cur=0.0f;for(int a=1;a<A;a++)if(best[xi][yj][zl][a]==k)cur=mems[xi][yj][zl][a];
                            ratio=(mems[xi][yj][zl][0]+ZERO)/(cur+ZERO);
                        }
                        float next=val+ratio;available[xi][yj][zl][k]=false;tree.addValue(next,xi,yj,zl,k);inserted.add(item(next,xi,yj,zl,k));
                    }
                }else available[x][y][z][k]=true;
            }
            events.add("["+item(val,x,y,z,k)+","+(regular==null?-1:(regular?1:0))+",["+String.join(",",inserted)+"]]");
        }
    }
    private static String labels(byte[][][] a,boolean skeleton){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(skeleton&&v==EMPTY?'a':(char)('0'+v));return b.toString();}
    private static String available(boolean[][][][] a){StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(a[x][y][z][k]?'1':'0');return b.toString();}
    private static String floats(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();}
    private static String ordering(float[][][] a){float[] flat=new float[N*N*N];int n=0;for(float[][] x:a)for(float[] y:x)for(float v:y)flat[n++]=v;return floats(flat);}
    private void emit(){
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.initGrowingLabels + propagateGrowing useApprox branch\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[9,9,9,2],\"events\":["+String.join(",",events)+"],\"segmentation\":\""+labels(seg,false)+"\",\"skeleton\":\""+labels(skel,true)+"\",\"available\":\""+available(available)+"\",\"ordering_bits\":"+ordering(ordering)+",\"volumes_bits\":"+floats(volume)+",\"queue_after\":[]}");
    }
    public static void main(String[] args){ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();ToadsApproximateGrowingPassReferenceHarness h=new ToadsApproximateGrowingPassReferenceHarness();h.initialize();h.propagate();h.emit();}
}
