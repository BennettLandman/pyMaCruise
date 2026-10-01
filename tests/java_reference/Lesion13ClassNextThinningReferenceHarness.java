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

/** Selected 13-class initThinningLabels + propagateThinning source excerpt.
 * The source topology-change branch is included and uses the selected 18/6 LUT.
 */
public final class Lesion13ClassNextThinningReferenceHarness {
    private static final int EMPTY = -1;
    private static final float ZERO = 1.0e-10f;
    private static final int[][] FACE = {{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private int nx,ny,nz,classes,lesionClass;
    private float volumeRatio,limit;
    private byte[] kinds;
    private byte[][][] seg,skel;
    private float[][][][] mem,shape;
    private float[][][] lesion,black,ordering;
    private boolean[][][][] available;
    private float[] volume,minimum;
    private BinaryHeap4D tree;
    private CriticalPointLUT lut;
    private final List<Item> initialQueue=new ArrayList<>();
    private final List<Event> events=new ArrayList<>();
    private static final class Item { final float v;final int x,y,z,k;Item(float v,int x,int y,int z,int k){this.v=v;this.x=x;this.y=y;this.z=z;this.k=k;} }
    private static final class Event { final Item pop;final boolean processed;final int candidate,topology;final List<Item> added;Event(Item p,boolean d,int c,int t,List<Item>a){pop=p;processed=d;candidate=c;topology=t;added=a;} }

    private void read(String path)throws Exception {
        try(DataInputStream in=new DataInputStream(new BufferedInputStream(new FileInputStream(path)))){
            nx=in.readInt();ny=in.readInt();nz=in.readInt();classes=in.readInt();lesionClass=in.readInt();volumeRatio=in.readFloat();limit=in.readFloat();
            kinds=new byte[classes];in.readFully(kinds);
            seg=new byte[nx][ny][nz];skel=new byte[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)seg[x][y][z]=in.readByte();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)skel[x][y][z]=in.readByte();
            mem=new float[nx][ny][nz][classes];shape=new float[nx][ny][nz][classes];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)mem[x][y][z][k]=in.readFloat();
            lesion=new float[nx][ny][nz];black=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)lesion[x][y][z]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)black[x][y][z]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)shape[x][y][z][k]=in.readFloat();
            available=new boolean[nx][ny][nz][classes];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)available[x][y][z][k]=in.readByte()!=0;
            ordering=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)ordering[x][y][z]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("trailing input");
        }
    }
    private float score(int x,int y,int z,int k){return kinds[k]==1||kinds[k]==2?mem[x][y][z][k]+lesion[x][y][z]+black[x][y][z]:mem[x][y][z][k];}
    private int best(int x,int y,int z){int id=0;float max=0.0f;for(int k=0;k<classes;k++){float value=score(x,y,z,k);if(value>max){id=k;max=mem[x][y][z][k];}}return id;}
    private float ratio(int x,int y,int z,int k,int b){return k==b?1.0f:(score(x,y,z,b)+ZERO)/(score(x,y,z,k)+ZERO);}
    private boolean boundary(int x,int y,int z,int k){for(int[] d:FACE)if(seg[x+d[0]][y+d[1]][z+d[2]]!=k)return true;return false;}
    private int pattern(boolean[][][] object){return lut.keyFromPattern(object,1,1,1);}
    private boolean[][][] neighborhood(int x,int y,int z,int a,int b,int c,boolean center){
        boolean[][][] object=new boolean[3][3][3];
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){
            int label=seg[x+i][y+j][z+l];
            object[1+i][1+j][1+l]=label==a||label==b||label==c;
        }
        object[1][1][1]=center;
        return object;
    }
    /** Direct excerpt of selected LesionToadSegmentation.isHomeomorphicSegmentation. */
    private boolean homeomorphic(int x,int y,int z,int k){
        int old=seg[x][y][z];
        if(old==k)return true;
        if(!lut.get(pattern(neighborhood(x,y,z,k,Integer.MIN_VALUE,Integer.MIN_VALUE,true))))return false;
        if(!lut.get(pattern(neighborhood(x,y,z,old,Integer.MIN_VALUE,Integer.MIN_VALUE,false))))return false;
        int[] labels=new int[26];int count=0;
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){
            if(i*i+j*j+l*l==0)continue;
            int label=seg[x+i][y+j][z+l];
            if(label==k||label==old)continue;
            boolean found=false;for(int n=0;n<count;n++)if(labels[n]==label){found=true;break;}
            if(!found)labels[count++]=label;
        }
        for(int n=0;n<count;n++){
            if(!lut.get(pattern(neighborhood(x,y,z,old,labels[n],Integer.MIN_VALUE,false))))return false;
        }
        for(int n=0;n<count;n++){
            if(!lut.get(pattern(neighborhood(x,y,z,k,labels[n],Integer.MIN_VALUE,true))))return false;
        }
        for(int n=0;n<count;n++)for(int m=n+1;m<count;m++){
            if(!lut.get(pattern(neighborhood(x,y,z,old,labels[n],labels[m],false))))return false;
        }
        for(int n=0;n<count;n++)for(int m=n+1;m<count;m++){
            if(!lut.get(pattern(neighborhood(x,y,z,k,labels[n],labels[m],true))))return false;
        }
        return true;
    }
    private void run(){
        volume=new float[classes];minimum=new float[classes];
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int current=seg[x][y][z];skel[x][y][z]=(byte)current;
            for(int k=0;k<classes;k++){available[x][y][z][k]=current==k;if(current==k)volume[k]+=k==lesionClass?mem[x][y][z][k]+lesion[x][y][z]+black[x][y][z]:mem[x][y][z][k];}
        }
        initialVolume=volume.clone();
        tree=new BinaryHeap4D(nx*ny*nz*classes,BinaryHeap4D.MAXTREE);
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){int k=seg[x][y][z];if(boundary(x,y,z,k)){available[x][y][z][k]=false;float p=-1.0f/ratio(x,y,z,k,best(x,y,z));tree.addValue(p,x,y,z,k);initialQueue.add(new Item(p,x,y,z,k));}}
        for(int k=0;k<classes;k++)minimum[k]=Numerics.max(100.0f,volumeRatio*volume[k]);
        float previous=0.0f;
        while(tree.isNotEmpty()&&previous>-limit){
            float v=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();tree.removeFirst();
            boolean did=false;int candidate=EMPTY,topology=-1;List<Item> inserted=new ArrayList<>();
            if(volume[k]>minimum[k]&&skel[x][y][z]==k){
                int id=best(x,y,z);boolean direct=false;int next=EMPTY;float bestPrior=0.0f;
                if(id!=k){
                    for(int[] d:FACE){int xi=x+d[0],yj=y+d[1],zl=z+d[2];
                        if(seg[xi][yj][zl]==id)direct=true;
                        else if(seg[xi][yj][zl]!=k&&shape[xi][yj][zl][id]>bestPrior){next=seg[xi][yj][zl];bestPrior=shape[xi][yj][zl][id];}
                    }
                    if(!direct&&next!=EMPTY)id=next;
                }
                boolean changeBoundary=id!=k&&(direct||next!=EMPTY);
                if(changeBoundary){candidate=id;changeBoundary=homeomorphic(x,y,z,id);topology=changeBoundary?1:0;}
                available[x][y][z][k]=false;ordering[x][y][z]=v;
                volume[k]-=k==lesionClass?mem[x][y][z][k]+lesion[x][y][z]+black[x][y][z]:mem[x][y][z][k];
                skel[x][y][z]=(byte)EMPTY;if(changeBoundary)seg[x][y][z]=(byte)id;did=true;
                for(short i=-1;i<=1;i++)for(short j=-1;j<=1;j++)for(short l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(available[xi][yj][zl][k]){float order=v-1.0f/ratio(xi,yj,zl,k,best(xi,yj,zl));if(volume[k]>minimum[k]&&order>-limit){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);inserted.add(new Item(order,xi,yj,zl,k));}}}
            }
            events.add(new Event(new Item(v,x,y,z,k),did,candidate,topology,inserted));previous=v;
        }
    }
    private static void item(DataOutputStream o,Item i)throws Exception{o.writeFloat(i.v);o.writeInt(i.x);o.writeInt(i.y);o.writeInt(i.z);o.writeInt(i.k);}
    private void write(String path)throws Exception{
        try(DataOutputStream o=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))){
            for(float v:initialVolume)o.writeFloat(v);for(float v:minimum)o.writeFloat(v);
            o.writeInt(initialQueue.size());for(Item i:initialQueue)item(o,i);
            o.writeInt(events.size());for(Event e:events){item(o,e.pop);o.writeByte(e.processed?1:0);o.writeInt(e.candidate);o.writeByte(e.topology);o.writeInt(e.added.size());for(Item i:e.added)item(o,i);}
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)o.writeByte(seg[x][y][z]);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)o.writeByte(skel[x][y][z]);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)o.writeByte(available[x][y][z][k]?1:0);
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)o.writeFloat(ordering[x][y][z]);
            for(float v:volume)o.writeFloat(v);
            List<Item> remain=new ArrayList<>();while(tree.isNotEmpty()){remain.add(new Item(tree.getFirst(),tree.getFirstX(),tree.getFirstY(),tree.getFirstZ(),tree.getFirstK()));tree.removeFirst();}
            o.writeInt(remain.size());for(Item i:remain)item(o,i);
        }
    }
    private float[] initialVolume;
    public static void main(String[] args)throws Exception{if(args.length!=2)throw new IllegalArgumentException("input.bin output.bin");Lesion13ClassNextThinningReferenceHarness h=new Lesion13ClassNextThinningReferenceHarness();h.read(args[0]);h.lut=new CriticalPointLUT("critical186LUT.raw.gz",200);if(!h.lut.loadCompressedPattern())throw new IllegalStateException("could not load 18/6 critical-point LUT");h.run();h.write(args[1]);}
}
