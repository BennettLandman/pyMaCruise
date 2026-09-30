import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Fixed-input excerpt for the selected continuing TOADS iteration:
 * initThinningLabels/propagateThinning queue followed by computeExactMemberships.
 * Queue transitions are deliberately restricted to the source branch where
 * every class volume is below minVolume=100; the harness asserts that premise.
 */
public final class ToadsNextDefaultMembershipReferenceHarness {
    private static final String TOADS_SHA="63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String HEAP_SHA="9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a";
    private static final int K=11, NC=2;
    private static final float ZERO=1.0e-10f, INF=1.0e10f, SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f;
    private static final float[][] CENTROID={{.28f,.31f,.34f,.37f,.40f,.43f,.46f,.49f,.52f,.55f,.58f},{.24f,.27f,.30f,.33f,.36f,.39f,.42f,.45f,.48f,.51f,.54f}};
    private static final float[][] VARIANCE={{.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},{1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}};
    private static final float[] MODALITY={.73f,.46f}, IMAGE_VARIANCE={.61f,.34f};
    private static final float[] EXP_CSF={.78f,.91f}, EXP_GM={1.16f,1.29f}, EXP_WM={1.55f,1.73f};
    private static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] SIM={
        {0,.11f,.14f,.18f,.13f,.16f,.12f,.15f,.09f,.10f,.08f},
        {.12f,0,.13f,.17f,.11f,.16f,.14f,.10f,.09f,.12f,.08f},
        {.15f,.12f,0,.13f,.18f,.12f,.16f,.11f,.10f,.07f,.09f},
        {.14f,.16f,.12f,0,.13f,.18f,.11f,.14f,.08f,.10f,.12f},
        {.11f,.13f,.17f,.12f,0,.11f,.19f,.14f,.10f,.09f,.07f},
        {.16f,.12f,.14f,.17f,.11f,0,.12f,.18f,.09f,.08f,.10f},
        {.13f,.15f,.11f,.14f,.18f,.12f,0,.13f,.08f,.09f,.11f},
        {.12f,.10f,.14f,.16f,.13f,.18f,.11f,0,.07f,.09f,.08f},
        {.10f,.09f,.12f,.08f,.11f,.09f,.10f,.07f,0,.13f,.16f},
        {.11f,.12f,.08f,.10f,.09f,.07f,.11f,.08f,.14f,0,.12f},
        {.08f,.10f,.09f,.12f,.07f,.11f,.10f,.09f,.17f,.13f,0}
    };
    private static int nx,ny,nz;
    private static float[][][][] image,field,bias,initial,shapes;
    private static short[][][][] labels;
    private static short[][][] segmentation;
    private static boolean[][] adjacency;
    private static short[][][] skeleton;
    private static boolean[][][][] available;
    private static PowerTable[][] powers;

    private static byte[] bytes(String path)throws Exception{ByteArrayOutputStream out=new ByteArrayOutputStream();try(GZIPInputStream in=new GZIPInputStream(new FileInputStream(path))){byte[] chunk=new byte[8192];int n;while((n=in.read(chunk))>=0)out.write(chunk,0,n);}return out.toByteArray();}
    private static ByteBuffer header(String path)throws Exception{return ByteBuffer.wrap(bytes(path)).order(ByteOrder.LITTLE_ENDIAN);}
    private static int dim(ByteBuffer b,int axis){return b.getShort(40+2*axis)&65535;}
    private static int offset(ByteBuffer b){return (int)b.getFloat(108);}
    private static float[][][][] readF4(String path,int channels)throws Exception{ByteBuffer b=header(path);if(b.getInt(0)!=348||b.getShort(70)!=16||b.getShort(72)!=32||b.getShort(40)!=4||dim(b,4)!=channels)throw new IllegalArgumentException("expected float32 NIfTI-1 four-volume: "+path);int x=dim(b,1),y=dim(b,2),z=dim(b,3);if(nx==0){nx=x;ny=y;nz=z;}if(x!=nx||y!=ny||z!=nz)throw new IllegalArgumentException("subject grid mismatch: "+path);float[][][][] out=new float[nx][ny][nz][channels];int o=offset(b);for(int c=0;c<channels;c++)for(int iz=0;iz<nz;iz++)for(int iy=0;iy<ny;iy++)for(int ix=0;ix<nx;ix++)out[ix][iy][iz][c]=b.getFloat(o+4*(((c*nz+iz)*ny+iy)*nx+ix));return out;}
    private static short[][][][] readS4(String path,int channels)throws Exception{ByteBuffer b=header(path);if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=4||dim(b,4)!=channels)throw new IllegalArgumentException("expected int16 NIfTI-1 four-volume: "+path);short[][][][] out=new short[nx][ny][nz][channels];int o=offset(b);for(int c=0;c<channels;c++)for(int iz=0;iz<nz;iz++)for(int iy=0;iy<ny;iy++)for(int ix=0;ix<nx;ix++)out[ix][iy][iz][c]=b.getShort(o+2*(((c*nz+iz)*ny+iy)*nx+ix));return out;}
    private static short[][][] readS3(String path)throws Exception{ByteBuffer b=header(path);if(b.getInt(0)!=348||b.getShort(70)!=4||b.getShort(72)!=16||b.getShort(40)!=3||dim(b,1)!=nx||dim(b,2)!=ny||dim(b,3)!=nz)throw new IllegalArgumentException("expected int16 NIfTI-1 segmentation");short[][][] out=new short[nx][ny][nz];int o=offset(b);for(int iz=0;iz<nz;iz++)for(int iy=0;iy<ny;iy++)for(int ix=0;ix<nx;ix++)out[ix][iy][iz]=b.getShort(o+2*((iz*ny+iy)*nx+ix));return out;}
    private static boolean direct(int x,int y,int z,int k){return labels[x][y][z][0]==k||labels[x][y][z][1]==k;}
    private static boolean neighbor(int x,int y,int z,int k){return labels[x][y][z][2]==k||labels[x][y][z][3]==k;}
    private static int category(int k){return k<=2?0:(k<=7?1:2);}
    private static float exponent(int t,int c){return t==0?EXP_CSF[c]:(t==1?EXP_GM[c]:EXP_WM[c]);}
    private static int[] bits(float[][][][] a){int[] out=new int[nx*ny*nz*a[0][0][0].length];int p=0;for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int c=0;c<a[x][y][z].length;c++)out[p++]=Float.floatToIntBits(a[x][y][z][c]);return out;}
    private static int[] bits(float[][][] a){int[] out=new int[nx*ny*nz];int p=0;for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)out[p++]=Float.floatToIntBits(a[x][y][z]);return out;}
    private static void append(StringBuilder s,int[] values){s.append('[');for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);}s.append(']');}
    private static float[][][][] dataEnergy(){float[][][][] out=new float[nx][ny][nz][K];powers=new PowerTable[3][NC];for(int t=0;t<3;t++)for(int c=0;c<NC;c++)powers[t][c]=new PowerTable(0.0f,4.0f,.000001f,exponent(t,c));for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++)for(int k=0;k<K;k++){boolean d=direct(x,y,z,k),n=neighbor(x,y,z,k);if(!d&&!n&&shapes[x][y][z][k]==0.0f)continue;float num=0.0f;if(k==0){for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num+=w*(image[x][y][z][c]==0.0f?0.0f:1.0f);}}else{int t=category(k);for(int c=0;c<NC;c++){float value=field[x][y][z][c]*image[x][y][z][c];value-=bias[x][y][z][c];float delta=value-CENTROID[c][k];float distance=delta*delta;float weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];num=(float)(num+(double)weight*powers[t][c].lookup((double)distance,(double)exponent(t,c)));}}out[x][y][z][k]=num;}return out;}
    private static float relationFactor(int x,int y,int z,int k){if(direct(x,y,z,k))return 1.0f;if(neighbor(x,y,z,k))return Math.max(.25f,Math.min(2.0f*FACTOR*shapes[x][y][z][k],1.0f));return Math.min(1.0f,2.0f*FACTOR*shapes[x][y][z][k]);}
    private static float[][][][] membership(float[][][][] energy,float[][][] denominator,float[] objective){float[][][][] out=new float[nx][ny][nz][K];for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<K;k++)out[x][y][z][k]=initial[x][y][z][k];float obj=0;int[][] ds={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){float den=0;float[] objectiveAt=new float[K];for(int k=0;k<K;k++){boolean d=direct(x,y,z,k),n=neighbor(x,y,z,k);float shape=shapes[x][y][z][k];if(!d&&!n&&shape==0){out[x][y][z][k]=0;continue;}float num=energy[x][y][z][k],e=num,neighborSum=0;for(int[] delta:ds)for(int m=0;m<K;m++)if(m!=k){float v=out[x+delta[0]][y+delta[1]][z+delta[2]][m];neighborSum+=v*v;}num+=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;e+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;float inverse=0;for(int m=0;m<K;m++)if(m!=k)inverse+=ATLAS_PRIOR[m]*SIM[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];num+=COEFFICIENT*inverse;e+=.5f*COEFFICIENT*inverse;float reciprocal=num>ZERO?1.0f/num:INF;out[x][y][z][k]=relationFactor(x,y,z,k)*reciprocal;den+=out[x][y][z][k];objectiveAt[k]=e;}denominator[x][y][z]=den;if(den>0)for(int k=0;k<K;k++){out[x][y][z][k]/=den;obj+=objectiveAt[k]*out[x][y][z][k]*out[x][y][z][k];}}objective[0]=obj;return out;}

    private static float bestRatio(int x,int y,int z,int k){int best=0;float max=0.0f;for(int m=0;m<K;m++)if(initial[x][y][z][m]>max){best=m;max=initial[x][y][z][m];}if(k==best)return 1.0f;return (initial[x][y][z][best]+ZERO)/(initial[x][y][z][k]+ZERO);}
    private static String queueResult(float volumeRatio,float limit){
        float[] volume=new float[K];for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){int k=segmentation[x][y][z];if(x>0&&x<nx-1&&y>0&&y<ny-1&&z>0&&z<nz-1)volume[k]+=initial[x][y][z][k];}
        float[] minimum=new float[K];for(int k=0;k<K;k++)minimum[k]=Math.max(100.0f,volumeRatio*volume[k]);
        skeleton=new short[nx][ny][nz];available=new boolean[nx][ny][nz][K];for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)skeleton[x][y][z]=-1;
        int[] count=new int[K];
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){int k=segmentation[x][y][z];skeleton[x][y][z]=(short)k;available[x][y][z][k]=true;}
        BinaryHeap4D queue=new BinaryHeap4D(nx*ny+ny*nz+nz*nx,BinaryHeap4D.MAXTREE);
        int[][] face={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){int k=segmentation[x][y][z];boolean boundary=false;for(int[] d:face)if(segmentation[x+d[0]][y+d[1]][z+d[2]]!=k)boundary=true;if(!boundary)continue;available[x][y][z][k]=false;boolean endpoint=false;int ni=x,nj=y,nl=z,neighborCount=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l>0){neighborCount++;ni=x+i;nj=y+j;nl=z+l;}if(neighborCount==1){neighborCount=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l>0&&segmentation[ni+i][nj+j][nl+l]==k)neighborCount++;endpoint=neighborCount<=2;}float order=endpoint?-0.0f:-1.0f/bestRatio(x,y,z,k);queue.addValue(order,x,y,z,k);count[k]++;}
        List<String> popped=new ArrayList<>();float value=0.0f;while(queue.isNotEmpty()&&value>-limit){int x=queue.getFirstX(),y=queue.getFirstY(),z=queue.getFirstZ(),k=queue.getFirstK();value=queue.getFirst();popped.add("["+Float.floatToIntBits(value)+","+x+","+y+","+z+","+k+"]");queue.removeFirst();if(!(volume[k]>minimum[k])||skeleton[x][y][z]!=k)continue;
            int target=0;float max=0;for(int m=0;m<K;m++)if(initial[x][y][z][m]>max){target=m;max=initial[x][y][z][m];}boolean direct=false;int next=-1;float bestPrior=-.1f;int candidate=target;if(target!=k)for(int[] d:face){int xx=x+d[0],yy=y+d[1],zz=z+d[2],neighborLabel=segmentation[xx][yy][zz];if(neighborLabel==target)direct=true;else if(neighborLabel!=k&&adjacency[neighborLabel][target]&&shapes[xx][yy][zz][target]>bestPrior){next=neighborLabel;bestPrior=shapes[xx][yy][zz][target];}}
            boolean change=target!=k&&(direct||next>=0);if(change&&!direct)candidate=next;if(change)change=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(segmentation,(short)x,(short)y,(short)z,(short)candidate);available[x][y][z][k]=false;volume[k]-=initial[x][y][z][k];skeleton[x][y][z]=-1;if(change)segmentation[x][y][z]=(short)candidate;
            for(int[] d:face){int xx=x+d[0],yy=y+d[1],zz=z+d[2];if(!available[xx][yy][zz][k])continue;int nnb=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(segmentation[xx+i][yy+j][zz+l]==k)nnb++;float order;if(nnb==1)order=value;else order=value-1.0f/bestRatio(xx,yy,zz,k);if(volume[k]>minimum[k]&&order>-limit){available[xx][yy][zz][k]=false;queue.addValue(order,xx,yy,zz,k);}}
        }
        StringBuilder s=new StringBuilder("{\"source\":\"MultipleToadSegmentation.initThinningLabels+propagateThinning\",\"source_sha256\":\"").append(TOADS_SHA).append("\",\"heap_source_sha256\":\"").append(HEAP_SHA).append("\",\"volumes_bits\":");append(s,toBits(volume));s.append(",\"minimum_bits\":");append(s,toBits(minimum));s.append(",\"boundary_counts\":");append(s,count);s.append(",\"popped\":[").append(String.join(",",popped)).append("],\"segmentation\":");append(s,bits(segmentation));s.append(",\"skeleton\":");append(s,bits(skeleton));s.append(",\"available\":");boolean first=true;s.append('[');for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<K;k++){if(!first)s.append(',');first=false;s.append(available[x][y][z][k]);}s.append(']');return s.append('}').toString();
    }
    private static int[] toBits(float[] a){int[] out=new int[a.length];for(int i=0;i<a.length;i++)out[i]=Float.floatToIntBits(a[i]);return out;}
    private static int[] bits(short[][][] a){int[] out=new int[nx*ny*nz];int p=0;for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)out[p++]=a[x][y][z];return out;}
    private static boolean[][] readAdj(String path)throws Exception{java.util.Scanner sc=new java.util.Scanner(new java.io.File(path));boolean[][] a=new boolean[K][K];int i=0;while(sc.hasNextLine()){String[] v=sc.nextLine().trim().split(",");for(int j=0;j<K;j++)a[i][j]=v[j].equals("1");i++;}sc.close();if(i!=K)throw new IllegalArgumentException("adjacency must be 11x11");return a;}
    public static void main(String[] args)throws Exception{if(args.length!=9)throw new IllegalArgumentException("args: images field psf memberships shapes closest segmentation adjacency volumeRatio");ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();image=readF4(args[0],NC);field=readF4(args[1],NC);bias=readF4(args[2],NC);initial=readF4(args[3],K);shapes=readF4(args[4],K);labels=readS4(args[5],4);segmentation=readS3(args[6]);adjacency=readAdj(args[7]);float volumeRatio=Float.parseFloat(args[8]),limit=1.0f; // selected fixture pins the wrapper's thinning limit
        float[][][][] energy=dataEnergy();float[][][] den=new float[nx][ny][nz];float[] obj={0};float[][][][] mem=membership(energy,den,obj);StringBuilder s=new StringBuilder("{\"algorithm_source_sha256\":\"").append(TOADS_SHA).append("\",\"queue\":").append(queueResult(volumeRatio,limit)).append(",\"class_energy_bits\":");append(s,bits(energy));s.append(",\"membership_bits\":");append(s,bits(mem));s.append(",\"denominator_bits\":");append(s,bits(den));s.append(",\"objective_bits\":").append(Float.floatToIntBits(obj[0])).append('}');System.out.println(s);}
}
