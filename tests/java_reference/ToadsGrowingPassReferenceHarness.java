import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Synthetic selected-source excerpt for initGrowingLabels + propagateGrowing. */
public final class ToadsGrowingPassReferenceHarness {
    private static final int N=9,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] seg=new byte[N][N][N],skel=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final float[] volume={1,2,3,4};
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MINTREE);
    private final List<String> events=new ArrayList<>();

    private void initialize(){
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            skel[x][y][z]=EMPTY;
            ordering[x][y][z]=-((x*13+y*7+z*3)%19)/10.0f;
            for(int k=0;k<C;k++)mems[x][y][z][k]=(((x*5+y*7+z*11+k*3)%29)+1)/30.0f;
            int center=2+((y+z)%3);
            if(y>=2&&y<=6&&z>=2&&z<=6&&Math.abs(x-center)<=1)seg[x][y][z]=1;
            if(x>=5&&x<=6&&y>=2&&y<=3&&z>=5&&z<=6)seg[x][y][z]=2;
            if(x==1&&y>=2&&y<=6&&z==7)seg[x][y][z]=3;
        }
        int c1=2+((4+3)%3);skel[c1][4][3]=1;
        skel[5][2][5]=2;
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int k=0;k<C;k++)available[x][y][z][k]=skel[x][y][z]==EMPTY;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++){
            if(!available[x][y][z][k])continue;
            boolean boundary=false;
            for(int[] d:FACE)if(skel[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            available[x][y][z][k]=false;
            int best=0;float max=0.0f;
            for(int c=0;c<C;c++)if(mems[x][y][z][c]>max){best=c;max=mems[x][y][z][c];}
            float ratio=k==best?1.0f:(mems[x][y][z][best]+ZERO)/(mems[x][y][z][k]+ZERO);
            tree.addValue(ordering[x][y][z]+ratio,x,y,z,k);
        }
    }
    private boolean homeomorphic(int x,int y,int z,int k){
        short[][][] a=new short[N][N][N];
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)a[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(a,(short)x,(short)y,(short)z,(short)k);
    }
    private static String item(float value,int x,int y,int z,int k){return "["+Float.floatToIntBits(value)+","+x+","+y+","+z+","+k+"]";}
    private void propagate(){
        while(tree.isNotEmpty()){
            float val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
            Boolean regular=null;List<String> inserted=new ArrayList<>();
            if(skel[x][y][z]==EMPTY){
                regular=homeomorphic(x,y,z,k);regular=homeomorphic(x,y,z,k);
                if(regular){
                    skel[x][y][z]=(byte)k;for(int m=0;m<C;m++)available[x][y][z][m]=false;
                    ordering[x][y][z]=val;volume[k]+=mems[x][y][z][k];seg[x][y][z]=(byte)k;
                    for(int[] d:FACE){int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;
                        int best=0;float max=0;for(int c=0;c<C;c++)if(mems[xi][yj][zl][c]>max){best=c;max=mems[xi][yj][zl][c];}
                        float ratio=k==best?1.0f:(mems[xi][yj][zl][best]+ZERO)/(mems[xi][yj][zl][k]+ZERO);
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
    private static String ints(short[][][] values){StringBuilder b=new StringBuilder("[");boolean first=true;for(short[][] x:values)for(short[] y:x)for(short v:y){if(!first)b.append(',');first=false;b.append(v);}return b.append(']').toString();}
    private static String ints(byte[][][][] values){StringBuilder b=new StringBuilder("[");boolean first=true;for(byte[][][] x:values)for(byte[][] y:x)for(byte[] z:y)for(byte v:z){if(!first)b.append(',');first=false;b.append(v);}return b.append(']').toString();}
    private static String bools(boolean[][] values){StringBuilder b=new StringBuilder("[");boolean first=true;for(boolean[] x:values)for(boolean v:x){if(!first)b.append(',');first=false;b.append(v?'1':'0');}return b.append(']').toString();}
    private String relationStateAfterGrowth(){
        int limit=2;
        ToadsRelationDiffusionReferenceHarness.Outcome relation=ToadsRelationDiffusionReferenceHarness.computeRelations(seg,C,mems,limit);
        byte[][][][] initialBoundary=ToadsBoundaryDistancesReferenceHarness.initial(N,N,N,C,true);
        ToadsBoundaryDistancesReferenceHarness.Outcome boundary=ToadsBoundaryDistancesReferenceHarness.compute(seg,initialBoundary,C,limit);
        byte[][][][] closest=ToadsClosestLabelsReferenceHarness.closestLabels(seg,limit);
        ToadsRelationListReferenceHarness.Result table=ToadsRelationListReferenceHarness.build(seg,C);
        return "{\"limit\":"+limit+",\"relation_map\":"+ints(relation.map)+",\"relation_membership\":"+bools(table.relations)
            +",\"relation_diffusion_iterations\":"+relation.iterations+",\"boundary_diffusion_passes\":"+boundary.passes
            +",\"initial_boundary_distances\":"+ints(initialBoundary)
            +",\"boundary_distances\":"+ints(boundary.distances)+",\"closest_labels\":"+ints(closest)+"}";
    }
    private void emit(){
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.initGrowingLabels + propagateGrowing queue exhaustion\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[9,9,9,4],\"events\":["+String.join(",",events)+"],\"segmentation\":\""+labels(seg,false)+"\",\"skeleton\":\""+labels(skel,true)+"\",\"available\":\""+available(available)+"\",\"ordering_bits\":"+ordering(ordering)+",\"volumes_bits\":"+floats(volume)+",\"queue_after\":[],\"relation_state_after_growth\":"+relationStateAfterGrowth()+"}");
    }
    public static void main(String[] args){ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();ToadsGrowingPassReferenceHarness h=new ToadsGrowingPassReferenceHarness();h.initialize();h.propagate();h.emit();}
}
