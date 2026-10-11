import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;

/**
 * Selected-source-order next-loop thinning slice after a registered TOADS
 * continuing branch. The formulas follow VandyCRUISE
 * MultipleToadSegmentation.initThinningLabels/propagateThinning (April 2006).
 */
public final class ToadsPostRegistrationNextThinningReferenceHarness {
    private static final float ZERO=1e-10f, LIMIT=3.0f, VOLUME_RATIO=.5f;
    private static final int EMPTY=-1;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final int nx,ny,nz,classes;
    private final byte[][][] seg,skel;
    private final float[][][][] mems,shape;
    private final boolean[][][][] available;
    private final boolean[][] relation;
    private final float[] volume,minimum,initialVolume;
    private final int[] boundaries;
    private final BinaryHeap4D tree;
    private final List<String> processed=new ArrayList<>();
    private final CriticalPointLUT lut=new CriticalPointLUT("/critical186LUT.raw.gz",200);

    private ToadsPostRegistrationNextThinningReferenceHarness(int x,int y,int z,int c) {
        nx=x;ny=y;nz=z;classes=c;
        seg=new byte[nx][ny][nz];skel=new byte[nx][ny][nz];
        mems=new float[nx][ny][nz][classes];shape=new float[nx][ny][nz][classes];
        available=new boolean[nx][ny][nz][classes];relation=new boolean[classes][classes];
        volume=new float[classes];minimum=new float[classes];initialVolume=new float[classes];
        boundaries=new int[classes];
        tree=new BinaryHeap4D(nx*ny*nz,BinaryHeap4D.MAXTREE);
    }

    private void read(String path) throws Exception {
        try(DataInputStream in=new DataInputStream(new FileInputStream(path))) {
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)seg[x][y][z]=in.readByte();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)mems[x][y][z][k]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)shape[x][y][z][k]=in.readFloat();
            for(int k=0;k<classes;k++)for(int j=0;j<classes;j++)relation[k][j]=in.readByte()!=0;
            // initThinningLabels overwrites the interior skeleton but keeps
            // the caller-owned border. Read the exact state carried by the
            // prior propagateGrowing call rather than inventing a new one.
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)skel[x][y][z]=in.readByte();
            if(in.read()!=-1)throw new IllegalArgumentException("unexpected trailing input bytes");
        }
    }

    private int bestMembership(int x,int y,int z) {
        int best=0;float max=0.0f;
        for(int k=0;k<classes;k++)if(mems[x][y][z][k]>max){best=k;max=mems[x][y][z][k];}
        return best;
    }
    private boolean homeomorphic(int x,int y,int z,int target) {
        short[][][] values=new short[nx][ny][nz];
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int l=0;l<nz;l++)values[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(values,(short)x,(short)y,(short)z,(short)target);
    }
    private void initThinningLabels() {
        for(int k=0;k<classes;k++)volume[k]=0.0f;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            int k=seg[x][y][z];skel[x][y][z]=(byte)k;
            for(int c=0;c<classes;c++) {
                available[x][y][z][c]=k==c;
                if(k==c)volume[c]+=mems[x][y][z][c];
            }
        }
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            int k=seg[x][y][z];boolean boundary=false;int neighbors=0,ni=0,nj=0,nl=0;
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++) {
                int d2=i*i+j*j+l*l;
                if(d2==1&&seg[x+i][y+j][z+l]!=k)boundary=true;
                else if(d2>0){neighbors++;ni=x+i;nj=y+j;nl=z+l;}
            }
            if(!boundary)continue;
            available[x][y][z][k]=false;
            boolean endpoint=false;
            if(neighbors==1) {
                neighbors=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(i*i+j*j+l*l>0&&seg[ni+i][nj+j][nl+l]==k)neighbors++;
                if(neighbors<=2)endpoint=true;
            }
            float order;
            if(endpoint)order=-0.0f;
            else {int best=bestMembership(x,y,z);float ratio=k==best?1.0f:(mems[x][y][z][best]+ZERO)/(mems[x][y][z][k]+ZERO);order=-1.0f/ratio;}
            tree.addValue(order,x,y,z,k);boundaries[k]++;
        }
        for(int k=0;k<classes;k++) {
            initialVolume[k]=volume[k];minimum[k]=Math.max(100.0f,VOLUME_RATIO*volume[k]);
        }
    }
    private void thin() {
        float val=0.0f;
        while(tree.isNotEmpty()&&val>-LIMIT) {
            val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();
            processed.add("["+Float.floatToIntBits(val)+","+x+","+y+","+z+","+k+"]");tree.removeFirst();
            if(volume[k]<=minimum[k]||skel[x][y][z]!=k)continue;
            int id=bestMembership(x,y,z),next=EMPTY;boolean direct=false;float bestPrior=-.1f;
            if(id!=k)for(int[] d:FACE) {
                int xi=x+d[0],yj=y+d[1],zl=z+d[2],neighbor=seg[xi][yj][zl];
                if(neighbor==id)direct=true;
                else if(neighbor!=k&&relation[neighbor][id]&&shape[xi][yj][zl][id]>bestPrior) {
                    next=neighbor;bestPrior=shape[xi][yj][zl][id];
                }
            }
            boolean change=id!=k&&(direct||next!=EMPTY);
            if(change&&!direct)id=next;
            if(change)change=homeomorphic(x,y,z,id);
            available[x][y][z][k]=false;volume[k]-=mems[x][y][z][k];skel[x][y][z]=EMPTY;
            if(change)seg[x][y][z]=(byte)id;
            for(int[] d:FACE) {
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
    private static String labels(byte[][][] values,boolean empty) {
        StringBuilder b=new StringBuilder();for(byte[][] x:values)for(byte[] y:x)for(byte v:y) {
            if(empty&&v==EMPTY)b.append('.');else b.append(v<10?(char)('0'+v):(char)('a'+v-10));
        }return b.toString();
    }
    private static String available(boolean[][][][] values) {
        StringBuilder b=new StringBuilder();for(boolean[][][] x:values)for(boolean[][] y:x)for(boolean[] z:y)for(boolean v:z)b.append(v?'1':'0');return b.toString();
    }
    private static String floats(float[] values) {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(values[i]));}return b.append(']').toString();
    }
    private static String ints(int[] values) {
        StringBuilder b=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)b.append(',');b.append(values[i]);}return b.append(']').toString();
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=5)throw new IllegalArgumentException("input.bin nx ny nz classes");
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        ToadsPostRegistrationNextThinningReferenceHarness h=new ToadsPostRegistrationNextThinningReferenceHarness(
                Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]));
        h.read(args[0]);h.initThinningLabels();h.thin();
        System.out.print("{\"source\":\"selected VandyCRUISE MultipleToadSegmentation.initThinningLabels + propagateThinning; short constructor useApprox=false\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"shape\":["+h.nx+","+h.ny+","+h.nz+","+h.classes+"],\"volume_ratio_bits\":"+Float.floatToIntBits(VOLUME_RATIO)+",\"limit_bits\":"+Float.floatToIntBits(LIMIT)+",\"initial_volumes_bits\":"+floats(h.initialVolume)+",\"minimum_volumes_bits\":"+floats(h.minimum)+",\"boundary_counts\":"+ints(h.boundaries)+",\"processed\":["+String.join(",",h.processed)+"],\"segmentation\":\""+labels(h.seg,false)+"\",\"skeleton\":\""+labels(h.skel,true)+"\",\"available\":\""+available(h.available)+"\",\"final_volumes_bits\":"+floats(h.volume)+"}");
    }
}
