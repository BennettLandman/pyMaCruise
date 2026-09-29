import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Selected-method formula excerpt for one LesionToadSegmentation growth event. */
public final class LesionGrowingEventReferenceHarness {
    private static final int N=5,C=4,EMPTY=-1;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final String kind;
    private final boolean useLesionMax;
    private final byte[][][] segmentation=new byte[N][N][N];
    private final byte[][][] skeleton=new byte[N][N][N];
    private final float[][][][] mems=new float[N][N][N][C];
    private final float[][][] lesionMembership=new float[N][N][N];
    private final float[][][] blackHoleMembership=new float[N][N][N];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final float[][][] ordering=new float[N][N][N];
    private final float[] volume={10.0f,20.0f,30.0f,40.0f};
    private final BinaryHeap4D tree=new BinaryHeap4D(32,BinaryHeap4D.MINTREE);
    private String[] objectTypes;
    private int lesionClass;
    private int eventLabel;
    private Boolean topologyAccepted;
    private java.util.List<int[]> inserted;

    private LesionGrowingEventReferenceHarness(String name) {
        kind=name;useLesionMax=name.equals("accepted_lesion_max");
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            segmentation[x][y][z]=0;skeleton[x][y][z]=0;
            for(int k=0;k<C;k++)mems[x][y][z][k]=((x*5+y*7+z*11+k*3)%29+1)/30.0f;
            lesionMembership[x][y][z]=((x+3*y+5*z)%7+1)/10.0f;
            blackHoleMembership[x][y][z]=((2*x+y+3*z)%9+1)/12.0f;
        }
        if(name.startsWith("accepted_")) {
            segmentation[1][2][2]=1;skeleton[1][2][2]=1;skeleton[2][2][2]=EMPTY;
            eventLabel=1;
            for(int[] d:FACE) {
                int x=2+d[0],y=2+d[1],z=2+d[2];
                if(skeleton[x][y][z]!=1){skeleton[x][y][z]=EMPTY;available[x][y][z][1]=true;}
            }
            lesionClass=name.equals("accepted_nonlesion")?0:1;
            objectTypes=name.equals("accepted_nonlesion")
                ?new String[]{"obj","obj","obj","obj"}
                :new String[]{"obj","lesion","blackhole","optimized_lesion"};
        } else if(name.equals("rejected")) {
            String encoded="020033121011213000023100323";int q=0;
            for(int x=1;x<4;x++)for(int y=1;y<4;y++)for(int z=1;z<4;z++)segmentation[x][y][z]=skeleton[x][y][z]=(byte)(encoded.charAt(q++)-'0');
            skeleton[2][2][2]=EMPTY;eventLabel=3;lesionClass=1;
            objectTypes=new String[]{"obj","lesion","blackhole","optimized_lesion"};
        } else {
            segmentation[2][2][2]=1;skeleton[2][2][2]=0;eventLabel=3;lesionClass=1;
            objectTypes=new String[]{"obj","lesion","blackhole","optimized_lesion"};
        }
        tree.addValue(0.25f,2,2,2,eventLabel);
        tree.addValue(0.5f,3,3,3,2);
    }

    private float bestAndRatio(int x,int y,int z,int k) {
        float value,maximum;
        int best=0;float bestMaximum=0.0f;
        for(int m=0;m<C;m++) {
            if(objectTypes[m].equals("lesion")||objectTypes[m].equals("blackhole")) {
                if(useLesionMax)value=edu.jhmi.rad.medic.utilities.Numerics.max(mems[x][y][z][m],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
                else value=mems[x][y][z][m]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
            } else value=mems[x][y][z][m];
            if(value>bestMaximum){best=m;bestMaximum=mems[x][y][z][m];}
        }
        if(k==best)value=1.0f;
        else {
            if(objectTypes[k].equals("lesion")||objectTypes[k].equals("blackhole")) {
                if(useLesionMax)value=edu.jhmi.rad.medic.utilities.Numerics.max(mems[x][y][z][k],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
                else value=mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
            } else value=mems[x][y][z][k];
            if(objectTypes[best].equals("lesion")||objectTypes[best].equals("blackhole")) {
                if(useLesionMax)maximum=edu.jhmi.rad.medic.utilities.Numerics.max(mems[x][y][z][best],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
                else maximum=mems[x][y][z][best]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
            } else maximum=mems[x][y][z][best];
            value=(maximum+1e-10f)/(value+1e-10f);
        }
        return value;
    }

    private void process() {
        int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();
        float val=tree.getFirst();tree.removeFirst();
        topologyAccepted=null;inserted=new java.util.ArrayList<>();
        if(skeleton[x][y][z]==EMPTY) {
            short[][][] local=new short[5][5][5];
            for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int l=0;l<5;l++)local[i][j][l]=segmentation[i][j][l];
            topologyAccepted=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(local,(short)2,(short)2,(short)2,(short)k);
            if(topologyAccepted) {
                skeleton[x][y][z]=(byte)k;
                for(int m=0;m<C;m++)available[x][y][z][m]=false;
                ordering[x][y][z]=val;
                if(k==lesionClass) {
                    if(useLesionMax)volume[k]+=edu.jhmi.rad.medic.utilities.Numerics.max(mems[x][y][z][k],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
                    else volume[k]+=mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
                } else volume[k]+=mems[x][y][z][k];
                segmentation[x][y][z]=(byte)k;
                for(int[] d:FACE) {
                    int xi=x+d[0],yj=y+d[1],zl=z+d[2];
                    if(available[xi][yj][zl][k]) {
                        float priority=val+bestAndRatio(xi,yj,zl,k);
                        available[xi][yj][zl][k]=false;
                        tree.addValue(priority,xi,yj,zl,k);
                        inserted.add(new int[]{Float.floatToIntBits(priority),xi,yj,zl,k});
                    }
                }
            } else available[x][y][z][k]=true;
        }
    }
    private static String labels(byte[][][] values,boolean skeleton) {
        StringBuilder b=new StringBuilder();for(byte[][] x:values)for(byte[] y:x)for(byte v:y)b.append(skeleton&&v==EMPTY?'a':(char)('0'+v));return b.toString();
    }
    private static String avail(boolean[][][][] values) {
        StringBuilder b=new StringBuilder();for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)b.append(values[x][y][z][k]?'1':'0');return b.toString();
    }
    private static String floats(float[] values) {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(values[i]));}return b.append(']').toString();
    }
    private String queue() {
        StringBuilder b=new StringBuilder("[");boolean first=true;
        while(tree.isNotEmpty()){if(!first)b.append(',');first=false;b.append('[').append(Float.floatToIntBits(tree.getFirst())).append(',').append(tree.getFirstX()).append(',').append(tree.getFirstY()).append(',').append(tree.getFirstZ()).append(',').append(tree.getFirstK()).append(']');tree.removeFirst();}
        return b.append(']').toString();
    }
    private String insertedJson() {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<inserted.size();i++){if(i>0)b.append(',');int[] d=inserted.get(i);b.append('[').append(d[0]).append(',').append(d[1]).append(',').append(d[2]).append(',').append(d[3]).append(',').append(d[4]).append(']');}return b.append(']').toString();
    }
    private void emit(boolean first) {
        if(!first)System.out.print(',');
        System.out.print("{\"name\":\""+kind+"\",\"label\":"+eventLabel+",\"lesion_class\":"+lesionClass+",\"use_lesion_max\":"+useLesionMax+",\"topology_accepted\":"+(topologyAccepted==null?"null":topologyAccepted.toString())+",\"segmentation\":\""+labels(segmentation,false)+"\",\"skeleton\":\""+labels(skeleton,true)+"\",\"available\":\""+avail(available)+"\",\"ordering_bits\":"+floats(flatOrdering())+",\"volume_bits\":"+floats(volume)+",\"inserted\":"+insertedJson()+",\"queue_after\":"+queue()+"}");
    }
    private float[] flatOrdering(){float[] out=new float[N*N*N];int i=0;for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)out[i++]=ordering[x][y][z];return out;}
    public static void main(String[] args) {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        String[] names={"accepted_lesion_sum","accepted_lesion_max","accepted_nonlesion","rejected","stale"};
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.propagateGrowing one item formula excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"numerics_sha256\":\"13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb\",\"cases\":[");
        for(int i=0;i<names.length;i++){LesionGrowingEventReferenceHarness h=new LesionGrowingEventReferenceHarness(names[i]);h.process();h.emit(i==0);}
        System.out.println("]}");
    }
}
