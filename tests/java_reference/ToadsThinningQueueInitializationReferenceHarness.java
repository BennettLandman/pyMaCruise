import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Source-formula excerpt for MultipleToadSegmentation.initThinningLabels. */
public final class ToadsThinningQueueInitializationReferenceHarness {
    private static final int NX=7,NY=7,NZ=7,C=4;
    private static final float ZERO=1e-10f;
    private static byte[][][] segmentation=new byte[NX][NY][NZ];
    private static byte[][][] skeleton=new byte[NX][NY][NZ];
    private static boolean[][][][] available=new boolean[NX][NY][NZ][C];
    private static float[][][][] mems=new float[NX][NY][NZ][C];
    private static float[] volume=new float[C];
    private static BinaryHeap4D tree=new BinaryHeap4D(NX*NY+NY*NZ+NZ*NX,BinaryHeap4D.MAXTREE);
    private static int[] nbound=new int[C];

    private static int bestMembership(int x,int y,int z){byte id=0;float max=0.0f;for(byte m=0;m<C;m++)if(mems[x][y][z][m]>max){id=m;max=mems[x][y][z][m];}return id;}
    private static float maxMembershipRatio(int x,int y,int z,int k,int b){if(k==b)return 1.0f;else return (mems[x][y][z][b]+ZERO)/(mems[x][y][z][k]+ZERO);}
    private static void initialize(){
        for(int k=0;k<C;k++)volume[k]=0.0f;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            segmentation[x][y][z]=(byte)((x+2*y+3*z)%C);
            skeleton[x][y][z]=(byte)-1;
            for(int k=0;k<C;k++){
                mems[x][y][z][k]=(((x*5+y*7+z*11+k*3)%29)+1)/30.0f;
                available[x][y][z][k]=(x==0||x==NX-1||y==0||y==NY-1||z==0||z==NZ-1)&&((x+y+z+k)%2==0);
            }
        }
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            skeleton[x][y][z]=segmentation[x][y][z];
            for(int k=0;k<C;k++){
                if(segmentation[x][y][z]==k){volume[k]+=mems[x][y][z][k];available[x][y][z][k]=true;}
                else available[x][y][z][k]=false;
            }
        }
        int[] lb=new int[6];float[] score=new float[6];int Nlb,best=0;boolean found;
        for(int k=0;k<C;k++)nbound[k]=0;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<C;k++){
            if(segmentation[x][y][z]==k){
                boolean isBoundary=false;Nlb=0;int Nngb=0;int ni=0,nj=0,nl=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){
                    if(i*i+j*j+l*l==1&&segmentation[x+i][y+j][z+l]!=k)isBoundary=true;
                    else if(i*i+j*j+l*l>0){Nngb++;ni=x+i;nj=y+j;nl=z+l;}
                }
                if(isBoundary){
                    available[x][y][z][k]=false;
                    boolean endpoint=false;
                    if(Nngb==1){Nngb=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l>0&&segmentation[ni+i][nj+j][nl+l]==k)Nngb++;if(Nngb<=2)endpoint=true;}
                    float order;
                    if(endpoint)order=-0.0f;
                    else order=-1.0f/maxMembershipRatio(x,y,z,k,bestMembership(x,y,z));
                    tree.addValue(order,x,y,z,k);nbound[k]++;
                }
            }
        }
    }
    private static String fbits(float v){return Integer.toString(Float.floatToIntBits(v));}
    private static void printQueue(){System.out.print('[');boolean first=true;while(tree.isNotEmpty()){if(!first)System.out.print(',');first=false;System.out.print('['+fbits(tree.getFirst())+','+tree.getFirstX()+','+tree.getFirstY()+','+tree.getFirstZ()+','+tree.getFirstK()+']');tree.removeFirst();}System.out.print(']');}
    private static void printVolume(){System.out.print('[');for(int k=0;k<C;k++){if(k>0)System.out.print(',');System.out.print(fbits(volume[k]));}System.out.print(']');}
    private static void printCounts(){System.out.print('[');for(int k=0;k<C;k++){if(k>0)System.out.print(',');System.out.print(nbound[k]);}System.out.print(']');}
    private static void printLabels(){System.out.print('[');boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){if(!first)System.out.print(',');first=false;System.out.print(skeleton[x][y][z]);}System.out.print(']');}
    private static void printAvailable(){System.out.print('[');boolean first=true;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<C;k++){if(!first)System.out.print(',');first=false;System.out.print(available[x][y][z][k]);}System.out.print(']');}
    public static void main(String[] args){
        initialize();
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.initThinningLabels\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"heap_source_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"shape\":[7,7,7,4],\"volumes_bits\":");printVolume();System.out.print(",\"boundary_counts\":");printCounts();System.out.print(",\"skeleton\":");printLabels();System.out.print(",\"available\":");printAvailable();System.out.print(",\"popped\":");printQueue();System.out.println('}');
    }
}
