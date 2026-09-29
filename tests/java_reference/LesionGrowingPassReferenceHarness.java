import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Selected LesionToadSegmentation growing-pass excerpt on a folded phantom. */
public final class LesionGrowingPassReferenceHarness {
    private static final int N=9,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final String name;
    private final boolean useLesionMax;
    private final byte[][][] seg=new byte[N][N][N],skel=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N],lesion=new float[N][N][N],blackHole=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final float[] volume={1,2,3,4};
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*N*C,BinaryHeap4D.MINTREE);
    private final List<String> events=new ArrayList<>();
    private final String[] objectTypes={"obj","lesion","blackhole","optimized_lesion"};

    private LesionGrowingPassReferenceHarness(String name) {
        this.name=name;
        this.useLesionMax=name.equals("max");
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            skel[x][y][z]=(byte)EMPTY;
            ordering[x][y][z]=-((x*13+y*7+z*3)%19)/10.0f;
            lesion[x][y][z]=((x+3*y+5*z)%7+1)/10.0f;
            blackHole[x][y][z]=((2*x+y+3*z)%9+1)/12.0f;
            for(int k=0;k<C;k++)mems[x][y][z][k]=((x*5+y*7+z*11+k*3)%29+1)/30.0f;
            int center=2+((y+z)%3);
            if(y>=2&&y<=6&&z>=2&&z<=6&&Math.abs(x-center)<=1)seg[x][y][z]=1;
            if(x>=5&&x<=6&&y>=2&&y<=3&&z>=5&&z<=6)seg[x][y][z]=2;
            if(x==1&&y>=2&&y<=6&&z==7)seg[x][y][z]=3;
        }
        int center=2+((4+3)%3);
        skel[center][4][3]=1;
        skel[5][2][5]=2;
        skel[1][4][7]=3;
    }

    private float score(int x,int y,int z,int k) {
        if(!objectTypes[k].equals("lesion")&&!objectTypes[k].equals("blackhole"))return mems[x][y][z][k];
        if(useLesionMax)return edu.jhmi.rad.medic.utilities.Numerics.max(mems[x][y][z][k],lesion[x][y][z],blackHole[x][y][z]);
        return mems[x][y][z][k]+lesion[x][y][z]+blackHole[x][y][z];
    }
    private int bestMembership(int x,int y,int z) {
        int best=0;float maximum=0.0f;
        for(int k=0;k<C;k++) {
            float candidate=score(x,y,z,k);
            if(candidate>maximum){best=k;maximum=mems[x][y][z][k];}
        }
        return best;
    }
    private float maxMembershipRatio(int x,int y,int z,int k,int best) {
        if(k==best)return 1.0f;
        return (score(x,y,z,best)+ZERO)/(score(x,y,z,k)+ZERO);
    }
    private void initialize() {
        tree.reset();tree.setMinTree();
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int k=0;k<C;k++)available[x][y][z][k]=skel[x][y][z]==EMPTY;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++) {
            if(!available[x][y][z][k])continue;
            boolean boundary=false;
            for(int[] d:FACE)if(skel[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            available[x][y][z][k]=false;
            int best=bestMembership(x,y,z);
            float ratio=maxMembershipRatio(x,y,z,k,best);
            tree.addValue(ordering[x][y][z]+ratio,x,y,z,k);
        }
    }
    private boolean homeomorphic(int x,int y,int z,int k) {
        short[][][] a=new short[N][N][N];
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)a[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(a,(short)x,(short)y,(short)z,(short)k);
    }
    private static String item(float value,int x,int y,int z,int k) {
        return "["+Float.floatToIntBits(value)+","+x+","+y+","+z+","+k+"]";
    }
    private void propagate() {
        while(tree.isNotEmpty()) {
            float value=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
            Boolean accepted=null;List<String> inserted=new ArrayList<>();
            if(skel[x][y][z]==EMPTY) {
                accepted=homeomorphic(x,y,z,k);
                if(accepted) {
                    skel[x][y][z]=(byte)k;
                    for(int m=0;m<C;m++)available[x][y][z][m]=false;
                    ordering[x][y][z]=value;
                    volume[k]+=k==1?score(x,y,z,k):mems[x][y][z][k];
                    seg[x][y][z]=(byte)k;
                    for(int[] d:FACE) {
                        int xi=x+d[0],yj=y+d[1],zl=z+d[2];
                        if(!available[xi][yj][zl][k])continue;
                        int best=bestMembership(xi,yj,zl);
                        float priority=value+maxMembershipRatio(xi,yj,zl,k,best);
                        available[xi][yj][zl][k]=false;
                        tree.addValue(priority,xi,yj,zl,k);
                        inserted.add(item(priority,xi,yj,zl,k));
                    }
                } else available[x][y][z][k]=true;
            }
            events.add("["+item(value,x,y,z,k)+","+(accepted==null?-1:(accepted?1:0))+",["+String.join(",",inserted)+"]]");
        }
    }
    private static String labels(byte[][][] a,boolean empty) {
        StringBuilder b=new StringBuilder();
        for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(empty&&v==EMPTY?'a':(char)('0'+v));
        return b.toString();
    }
    private static String available(boolean[][][][] a) {
        StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(a[x][y][z][k]?'1':'0');return b.toString();
    }
    private static String floats(float[] a) {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();
    }
    private static String ordering(float[][][] a) {
        float[] out=new float[N*N*N];int i=0;for(float[][] x:a)for(float[] y:x)for(float v:y)out[i++]=v;return floats(out);
    }
    private void emit(boolean first) {
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+name+"\",\"use_lesion_max\":"+useLesionMax+",\"initial_segmentation\":\""+labels(initialSegmentation(),false)+"\",\"initial_skeleton\":\""+labels(initialSkeleton(),true)+"\",\"events\":["+String.join(",",events)+"],\"segmentation\":\""+labels(seg,false)+"\",\"skeleton\":\""+labels(skel,true)+"\",\"available\":\""+available(available)+"\",\"ordering_bits\":"+ordering(ordering)+",\"volume_bits\":"+floats(volume)+",\"queue_after\":[]}");
    }
    private byte[][][] initialSegmentation() {byte[][][] a=new byte[N][N][N];for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){int center=2+((y+z)%3);if(y>=2&&y<=6&&z>=2&&z<=6&&Math.abs(x-center)<=1)a[x][y][z]=1;if(x>=5&&x<=6&&y>=2&&y<=3&&z>=5&&z<=6)a[x][y][z]=2;if(x==1&&y>=2&&y<=6&&z==7)a[x][y][z]=3;}return a;}
    private byte[][][] initialSkeleton() {byte[][][] a=new byte[N][N][N];for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)a[x][y][z]=(byte)EMPTY;int center=2+((4+3)%3);a[center][4][3]=1;a[5][2][5]=2;a[1][4][7]=3;return a;}
    public static void main(String[] args) {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.initGrowingLabels + propagateGrowing queue exhaustion, selected-method formula harness\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"cases\":[");
        String[] names={"sum","max"};
        for(int i=0;i<names.length;i++){LesionGrowingPassReferenceHarness h=new LesionGrowingPassReferenceHarness(names[i]);h.initialize();h.propagate();h.emit(i==0);}
        System.out.println("]}");
    }
}
