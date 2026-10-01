import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;
import edu.jhmi.rad.medic.utilities.Numerics;

/** Selected LesionToadSegmentation growing methods on caller-supplied state. */
public final class LesionGrowingPassDynamicReferenceHarness {
    private static final int EMPTY = -1;
    private static final float ZERO = 1.0e-10f;
    private static final int[][] FACE = {{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final CriticalPointLUT LUT = new CriticalPointLUT("/edu/jhmi/rad/medic/structures/critical186LUT.raw.gz", 200);

    private int nx, ny, nz, classes, lesionClass;
    private boolean useLesionMax;
    private byte[] objectKinds; // 1=LESION, 2=BLACKHOLE, otherwise ordinary
    private byte[][][] segmentation, skeleton;
    private float[][][][] memberships;
    private boolean[][][][] available;
    private float[][][] ordering, lesionMembership, blackHoleMembership;
    private float[] volume;
    private BinaryHeap4D tree;
    private final List<Event> events = new ArrayList<>();

    private static final class Item {
        final float value; final int x,y,z,k;
        Item(float value,int x,int y,int z,int k){this.value=value;this.x=x;this.y=y;this.z=z;this.k=k;}
    }
    private static final class Event {
        final Item popped; final int accepted; final List<Item> inserted;
        Event(Item popped,int accepted,List<Item> inserted){this.popped=popped;this.accepted=accepted;this.inserted=inserted;}
    }

    private void read(String path) throws Exception {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path)))) {
            nx=in.readInt(); ny=in.readInt(); nz=in.readInt(); classes=in.readInt(); lesionClass=in.readInt(); useLesionMax=in.readInt()!=0;
            objectKinds=new byte[classes]; in.readFully(objectKinds);
            segmentation=new byte[nx][ny][nz]; skeleton=new byte[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)segmentation[x][y][z]=in.readByte();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)skeleton[x][y][z]=in.readByte();
            memberships=new float[nx][ny][nz][classes];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)memberships[x][y][z][k]=in.readFloat();
            ordering=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)ordering[x][y][z]=in.readFloat();
            volume=new float[classes]; for(int k=0;k<classes;k++)volume[k]=in.readFloat();
            lesionMembership=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)lesionMembership[x][y][z]=in.readFloat();
            blackHoleMembership=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)blackHoleMembership[x][y][z]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("unexpected trailing input");
        }
    }

    private float score(int x,int y,int z,int k) {
        if(objectKinds[k]!=1 && objectKinds[k]!=2)return memberships[x][y][z][k];
        if(useLesionMax)return Numerics.max(memberships[x][y][z][k],lesionMembership[x][y][z],blackHoleMembership[x][y][z]);
        return memberships[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z];
    }
    private int bestMembership(int x,int y,int z) {
        int best=0;float maximum=0.0f;
        for(int k=0;k<classes;k++) { float candidate=score(x,y,z,k); if(candidate>maximum){best=k;maximum=memberships[x][y][z][k];} }
        return best;
    }
    private float maxMembershipRatio(int x,int y,int z,int k,int best) {
        if(k==best)return 1.0f;
        return (score(x,y,z,best)+ZERO)/(score(x,y,z,k)+ZERO);
    }
    private void initialize() {
        available=new boolean[nx][ny][nz][classes];
        tree=new BinaryHeap4D(nx*ny*nz*classes,BinaryHeap4D.MINTREE); tree.reset();tree.setMinTree();
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++)for(int k=0;k<classes;k++)available[x][y][z][k]=skeleton[x][y][z]==EMPTY;
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++) {
            if(!available[x][y][z][k])continue;
            boolean boundary=false;for(int[] d:FACE)if(skeleton[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;available[x][y][z][k]=false;
            tree.addValue(ordering[x][y][z]+maxMembershipRatio(x,y,z,k,bestMembership(x,y,z)),x,y,z,k);
        }
    }
    private boolean homeomorphic(int x,int y,int z,int k) {
        if(segmentation[x][y][z]==k)return true;
        boolean[][][] obj=new boolean[3][3][3];
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)obj[1+i][1+j][1+l]=segmentation[x+i][y+j][z+l]==k;
        obj[1][1][1]=true;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)obj[1+i][1+j][1+l]=segmentation[x+i][y+j][z+l]==segmentation[x][y][z];
        obj[1][1][1]=false;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        int count=0;short[] labels=new short[26];
        for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l>0&&segmentation[x+i][y+j][z+l]!=k&&segmentation[x+i][y+j][z+l]!=segmentation[x][y][z]){
            boolean found=false;for(int n=0;n<count;n++)if(segmentation[x+i][y+j][z+l]==labels[n]){found=true;break;}if(!found)labels[count++]=segmentation[x+i][y+j][z+l];
        }
        for(int n=0;n<count;n++){
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)obj[1+i][1+j][1+l]=segmentation[x+i][y+j][z+l]==segmentation[x][y][z]||segmentation[x+i][y+j][z+l]==labels[n];
            obj[1][1][1]=false;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        }
        for(int n=0;n<count;n++){
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)obj[1+i][1+j][1+l]=segmentation[x+i][y+j][z+l]==k||segmentation[x+i][y+j][z+l]==labels[n];
            obj[1][1][1]=true;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        }
        for(int n=0;n<count;n++)for(int m=n+1;m<count;m++){
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){short v=segmentation[x+i][y+j][z+l];obj[1+i][1+j][1+l]=v==segmentation[x][y][z]||v==labels[n]||v==labels[m];}
            obj[1][1][1]=false;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        }
        for(int n=0;n<count;n++)for(int m=n+1;m<count;m++){
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){short v=segmentation[x+i][y+j][z+l];obj[1+i][1+j][1+l]=v==k||v==labels[n]||v==labels[m];}
            obj[1][1][1]=true;if(!LUT.get(LUT.keyFromPattern(obj,1,1,1)))return false;
        }
        return true;
    }
    private void propagate() {
        while(tree.isNotEmpty()) {
            float value=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
            int accepted=-1;List<Item> inserted=new ArrayList<>();
            if(skeleton[x][y][z]==EMPTY) {
                boolean regular=homeomorphic(x,y,z,k);accepted=regular?1:0;
                if(regular){
                    skeleton[x][y][z]=(byte)k;for(int m=0;m<classes;m++)available[x][y][z][m]=false;ordering[x][y][z]=value;
                    float increment=k==lesionClass?(useLesionMax?Numerics.max(memberships[x][y][z][k],lesionMembership[x][y][z],blackHoleMembership[x][y][z]):memberships[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMembership[x][y][z]):memberships[x][y][z][k];
                    volume[k]+=increment;segmentation[x][y][z]=(byte)k;
                    for(int[] d:FACE){int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;int best=bestMembership(xi,yj,zl);float priority=value+maxMembershipRatio(xi,yj,zl,k,best);available[xi][yj][zl][k]=false;tree.addValue(priority,xi,yj,zl,k);inserted.add(new Item(priority,xi,yj,zl,k));}
                } else available[x][y][z][k]=true;
            }
            events.add(new Event(new Item(value,x,y,z,k),accepted,inserted));
        }
    }
    private void write(String path) throws Exception {
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))){
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)out.writeByte(segmentation[x][y][z]);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)out.writeByte(skeleton[x][y][z]);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)out.writeByte(available[x][y][z][k]?1:0);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)out.writeFloat(ordering[x][y][z]);
            for(float v:volume)out.writeFloat(v);
            out.writeInt(events.size());
            for(Event e:events){writeItem(out,e.popped);out.writeByte(e.accepted);out.writeInt(e.inserted.size());for(Item i:e.inserted)writeItem(out,i);}
            if(out instanceof DataOutputStream){} // keep the output format explicit and binary
        }
    }
    private static void writeItem(DataOutputStream out,Item i)throws Exception{out.writeFloat(i.value);out.writeInt(i.x);out.writeInt(i.y);out.writeInt(i.z);out.writeInt(i.k);}
    public static void main(String[] args)throws Exception {
        if(args.length!=2)throw new IllegalArgumentException("input.bin output.bin required");
        if(!LUT.loadCompressedPattern())throw new IllegalStateException("could not load 18/6 LUT");
        LesionGrowingPassDynamicReferenceHarness h=new LesionGrowingPassDynamicReferenceHarness();h.read(args[0]);h.initialize();h.propagate();h.write(args[1]);
    }
}
