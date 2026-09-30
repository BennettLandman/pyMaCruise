import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** One selected LesionTOADS heap pop composed with local removal and neighbor scheduling. */
public final class LesionToadsComposedThinningEventReferenceHarness {
    private static final int N=5,C=3,X=2,Y=2,Z=2,K=0;
    private static final String NEIGHBORHOOD="212221122120200120112112212";
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final short[][][] seg=new short[N][N][N], skeleton=new short[N][N][N];
    private static final float[][][][] mem=new float[N][N][N][C], shape=new float[N][N][N][C];
    private static final float[][][] lesion=new float[N][N][N], black=new float[N][N][N], ordering=new float[N][N][N];
    private static final boolean[][][][] available=new boolean[N][N][N][C];
    private static final float[] volume={12.0f,17.0f,9.0f}, minimum={2.0f,2.0f,2.0f};

    private static float score(int x,int y,int z,int c){return c==1?(mem[x][y][z][c]+lesion[x][y][z])+black[x][y][z]:mem[x][y][z][c];}
    private static int best(int x,int y,int z){int id=0;float max=0.0f;for(int c=0;c<C;c++){float s=score(x,y,z,c);if(s>max){id=c;max=mem[x][y][z][c];}}return id;}
    private static float ratio(int x,int y,int z,int k,int b){if(k==b)return 1.0f;return (score(x,y,z,b)+1e-10f)/(score(x,y,z,k)+1e-10f);}
    private static String bits(float value){return Integer.toString(Float.floatToRawIntBits(value));}
    private static void append(StringBuilder b,String v){if(b.length()>0)b.append(',');b.append(v);}
    public static void main(String[] args){
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){
            seg[x][y][z]=0;skeleton[x][y][z]=0;
            for(int c=0;c<C;c++){mem[x][y][z][c]=c==0?0.125f:(c==1?0.25f:0.0625f);shape[x][y][z][c]=0.0f;}
            lesion[x][y][z]=0.5f;black[x][y][z]=0.25f;ordering[x][y][z]=0.0f;
        }
        for(int n=0;n<27;n++){int x=n/9+1,y=(n/3)%3+1,z=n%3+1;seg[x][y][z]=(short)(NEIGHBORHOOD.charAt(n)-'0');skeleton[x][y][z]=seg[x][y][z];}
        skeleton[X][Y][Z]=K;ordering[X][Y][Z]=0.5f;
        for(int[] d:FACE)available[X+d[0]][Y+d[1]][Z+d[2]][K]=true;

        BinaryHeap4D tree=new BinaryHeap4D(16,BinaryHeap4D.MAXTREE);
        tree.addValue(-0.25f,X,Y,Z,K);tree.addValue(-0.5f,4,4,4,2);
        List<String> inserted=new ArrayList<>();
        float val=tree.getFirst();short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
        int id=best(x,y,z),candidate=-1;boolean changed=false,topology=false;
        if(id!=k){boolean direct=false;float bestPrior=0.0f;int next=-1;
            for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(seg[xi][yj][zl]==id)direct=true;else if(seg[xi][yj][zl]!=k&&shape[xi][yj][zl][id]>bestPrior){next=seg[xi][yj][zl];bestPrior=shape[xi][yj][zl][id];}}
            if(direct){candidate=id;changed=true;}else if(next!=-1){candidate=next;changed=true;}
            if(changed){topology=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(seg,x,y,z,(short)candidate);changed=topology;}
        }
        if(volume[k]>minimum[k]){
            available[x][y][z][k]=false;ordering[x][y][z]=val;volume[k]-=mem[x][y][z][k];skeleton[x][y][z]=-1;if(changed)seg[x][y][z]=(short)candidate;
            for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(available[xi][yj][zl][k]){float order=val-1.0f/ratio(xi,yj,zl,k,best(xi,yj,zl));if(volume[k]>minimum[k]&&order>-100.0f){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);inserted.add("["+bits(order)+","+xi+","+yj+","+zl+","+k+"]");}}}
        }
        StringBuilder ins=new StringBuilder();for(String item:inserted)append(ins,item);
        System.out.print("{\"popped\":["+bits(val)+","+x+","+y+","+z+","+k+"],\"best\":"+id+",\"candidate\":"+candidate+",\"topology\":"+topology+",\"seg\":"+seg[x][y][z]+",\"skeleton\":"+skeleton[x][y][z]+",\"available_center\":["+available[x][y][z][0]+","+available[x][y][z][1]+","+available[x][y][z][2]+"],\"ordering_bits\":"+bits(ordering[x][y][z])+",\"volume_bits\":["+bits(volume[0])+","+bits(volume[1])+","+bits(volume[2]) +"],\"inserted\":["+ins+"],\"remaining\":[");
        boolean first=true;while(tree.isNotEmpty()){if(!first)System.out.print(',');first=false;System.out.print("["+bits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");tree.removeFirst();}System.out.println("]}");
    }
}
