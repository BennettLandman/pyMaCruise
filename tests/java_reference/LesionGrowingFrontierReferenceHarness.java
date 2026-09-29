import java.util.ArrayList;
import java.util.List;

import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected-method formula excerpt for LesionToadSegmentation growing setup. */
public final class LesionGrowingFrontierReferenceHarness {
    private static final int N=7,C=4,EMPTY=-1;
    private static final float ZERO=1e-10f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] skeleton=new byte[N][N][N];
    private final float[][][] ordering=new float[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final float[][][] lesionMembership=new float[N][N][N];
    private final float[][][] blackHoleMembership=new float[N][N][N];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MINTREE);
    private final String kind;
    private final boolean useLesionMax;
    private final String[] objectTypes;

    private LesionGrowingFrontierReferenceHarness(String name) {
        kind=name;
        useLesionMax=name.equals("max");
        objectTypes=name.equals("optimized_types")
            ? new String[]{"optimized_lesion","optimized_blackHole","obj","mask"}
            : new String[]{"obj","lesion","blackhole","optimized_lesion"};
    }
    private void initializeInputs() {
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            skeleton[x][y][z]=EMPTY;
            ordering[x][y][z]=-((x*13+y*7+z*3)%19)/10.0f;
            for(int k=0;k<C;k++) {
                int numerator=(x*5+y*7+z*11+k*3)%29+1;
                mems[x][y][z][k]=kind.equals("zero_ties")?0.0f:numerator/30.0f;
            }
            lesionMembership[x][y][z]=kind.equals("zero_ties")?0.0f:((x+3*y+5*z)%7+1)/10.0f;
            blackHoleMembership[x][y][z]=kind.equals("zero_ties")?0.0f:((2*x+y+3*z)%9+1)/12.0f;
        }
        skeleton[3][3][2]=0;skeleton[3][3][4]=1;
        skeleton[2][3][3]=2;skeleton[4][3][3]=3;
    }
    // Exact source behavior from bestMembership, including the running-max assignment.
    private int bestMembership(int x,int y,int z) {
        int id=0;float max=0.0f;
        for(int m=0;m<C;m++) {
            float mem;
            if(objectTypes[m].equals("lesion")||objectTypes[m].equals("blackhole")) {
                if(useLesionMax)mem=Numerics.max(mems[x][y][z][m],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
                else mem=mems[x][y][z][m]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
            } else mem=mems[x][y][z][m];
            if(mem>max) {
                id=m;
                max=mems[x][y][z][m]; // source quirk: stores ordinary membership
            }
        }
        return id;
    }
    private float maxMembershipRatio(int x,int y,int z,int k,int best) {
        if(k==best)return 1.0f;
        float cur,max;
        if(objectTypes[k].equals("lesion")||objectTypes[k].equals("blackhole")) {
            if(useLesionMax)cur=Numerics.max(mems[x][y][z][k],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
            else cur=mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
        } else cur=mems[x][y][z][k];
        if(objectTypes[best].equals("lesion")||objectTypes[best].equals("blackhole")) {
            if(useLesionMax)max=Numerics.max(mems[x][y][z][best],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
            else max=mems[x][y][z][best]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
        } else max=mems[x][y][z][best];
        return (max+ZERO)/(cur+ZERO);
    }
    private void initialize() {
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int n=0;n<C;n++)available[x][y][z][n]=skeleton[x][y][z]==EMPTY;
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(byte k=0;k<C;k++) {
            if(!available[x][y][z][k])continue;
            boolean boundary=false;
            for(int[] d:FACE)if(skeleton[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            available[x][y][z][k]=false;
            int best=bestMembership(x,y,z);
            float order=ordering[x][y][z]+maxMembershipRatio(x,y,z,k,best);
            tree.addValue(order,x,y,z,k);
        }
    }
    private static String availableBits(boolean[][][][] values) {
        StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(values[x][y][z][k]?'1':'0');return b.toString();
    }
    private String queue() {
        List<String> out=new ArrayList<>();
        while(tree.isNotEmpty()) {
            out.add("["+Float.floatToIntBits(tree.getFirst())+","+tree.getFirstX()+","+tree.getFirstY()+","+tree.getFirstZ()+","+tree.getFirstK()+"]");
            tree.removeFirst();
        }
        return "["+String.join(",",out)+"]";
    }
    private void emit(boolean first) {
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+kind+"\",\"use_lesion_max\":"+useLesionMax+",\"object_types\":[\""+String.join("\",\"",objectTypes)+"\"],\"available\":\""+availableBits(available)+"\",\"queue\":"+queue()+"}");
    }
    public static void main(String[] args) {
        String[] cases={"additive","max","zero_ties","optimized_types"};
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.initGrowingLabels + bestMembership + maxMembershipRatio formula excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"numerics_sha256\":\"13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb\",\"shape\":[7,7,7,4],\"zero\":1.0e-10,\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"cases\":[");
        for(int i=0;i<cases.length;i++){LesionGrowingFrontierReferenceHarness h=new LesionGrowingFrontierReferenceHarness(cases[i]);h.initializeInputs();h.initialize();h.emit(i==0);}
        System.out.println("]}");
    }
}
