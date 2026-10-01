import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;

import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/**
 * Selected-method composition for the first exact-membership TOADS thinning
 * pass on the MedicAlgorithmToads08 default constructor path. The formulas
 * follow MultipleToadPreprocess.initialExactMemberships and
 * MultipleToadSegmentation.initThinningLabels/propagateThinning.
 */
public final class ToadsDefaultFirstThinningReferenceHarness {
    private static final int N=36, C=11, BORDER=10, EMPTY=-1;
    private static final float ZERO=1.0e-20f, THIN_ZERO=1.0e-10f;
    private static final float LIMIT=3.0f, VOLUME_RATIO=0.5f;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private final byte[][][] seg=new byte[N][N][N], skel=new byte[N][N][N];
    private final float[][][][] image=new float[2][N][N][N];
    private final float[][] centroid=new float[2][C];
    private final float[][][][] mems=new float[N][N][N][C];
    private final float[][][][] shape=new float[N][N][N][C];
    private final boolean[][][][] available=new boolean[N][N][N][C];
    private final boolean[][] relation=new boolean[C][C];
    private final float[] volume=new float[C], initialVolume=new float[C], minimum=new float[C];
    private final int[] boundaries=new int[C];
    private final BinaryHeap4D tree=new BinaryHeap4D(N*N*3,BinaryHeap4D.MAXTREE);
    private final List<String> processed=new ArrayList<>();

    private static int label(int x,int y,int z) {
        return Math.floorMod(x+2*y+3*z+(y*z)%3+(x*y)%2,C);
    }
    private static float intensity(int channel,int x,int y,int z,int k) {
        if(channel==0) {
            int delta=(x*13+y*7+z*3)%17-8;
            return (k+1)/12.0f+delta/1700.0f;
        }
        int delta=(x*5+y*11+z*7)%19-9;
        return (11-k)/12.0f+delta/1900.0f;
    }
    private void initializeMembershipsAndState() {
        for(int k=0;k<C;k++) {centroid[0][k]=(k+1)/12.0f;centroid[1][k]=(11-k)/12.0f;}
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            int current=label(x,y,z); seg[x][y][z]=(byte)current;skel[x][y][z]=(byte)EMPTY;
            image[0][x][y][z]=intensity(0,x,y,z,current);
            image[1][x][y][z]=intensity(1,x,y,z,current);
        }
        // MultipleToadPreprocess.initialExactMemberships. On this cropped
        // synthetic image, x0i/y0i/z0i are zero and bs is the source border.
        for(int x=BORDER;x<N-BORDER-1;x++)for(int y=BORDER;y<N-BORDER-1;y++)for(int z=BORDER;z<N-BORDER-1;z++) {
            float den=0.0f;
            for(int k=0;k<C;k++) {
                float num=0.0f;
                for(int c=0;c<2;c++) {
                    float difference=image[c][x][y][z]-centroid[c][k];
                    num+=difference*difference;
                }
                if(num>ZERO)num=1.0f/num;else num=1.0e20f;
                mems[x][y][z][k]=num;den+=num;
            }
            if(den>0.0f)for(int k=0;k<C;k++)mems[x][y][z][k]=mems[x][y][z][k]/den;
        }

        // initThinningLabels: reset volumes, copy labels into the skeleton,
        // initialize per-class availability, then seed boundary candidates.
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int k=seg[x][y][z];skel[x][y][z]=(byte)k;volume[k]+=mems[x][y][z][k];
            for(int c=0;c<C;c++)available[x][y][z][c]=c==k;
        }
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int k=seg[x][y][z];boolean boundary=false;int neighborCount=0,ni=0,nj=0,nl=0;
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++) {
                int d2=i*i+j*j+l*l;
                if(d2==1&&seg[x+i][y+j][z+l]!=k)boundary=true;
                else if(d2>0){neighborCount++;ni=x+i;nj=y+j;nl=z+l;}
            }
            if(!boundary)continue;
            available[x][y][z][k]=false;
            boolean endpoint=false;
            if(neighborCount==1) {
                neighborCount=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(i*i+j*j+l*l>0&&seg[ni+i][nj+j][nl+l]==k)neighborCount++;
                endpoint=neighborCount<=2;
            }
            float order;
            if(endpoint)order=-0.0f;
            else {
                int best=0;float max=0.0f;
                for(int c=0;c<C;c++)if(mems[x][y][z][c]>max){best=c;max=mems[x][y][z][c];}
                float ratio=k==best?1.0f:(mems[x][y][z][best]+THIN_ZERO)/(mems[x][y][z][k]+THIN_ZERO);
                order=-1.0f/ratio;
            }
            tree.addValue(order,x,y,z,k);boundaries[k]++;
        }
        System.arraycopy(volume,0,initialVolume,0,C);
        for(int k=0;k<C;k++)minimum[k]=Math.max(100.0f,VOLUME_RATIO*volume[k]);
    }
    private int bestMembership(int x,int y,int z) {
        int best=0;float max=0.0f;
        for(int k=0;k<C;k++)if(mems[x][y][z][k]>max){best=k;max=mems[x][y][z][k];}
        return best;
    }
    private boolean topology(int x,int y,int z,int target) {
        short[][][] copy=new short[N][N][N];
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)copy[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(copy,(short)x,(short)y,(short)z,(short)target);
    }
    private void oneEvent() {
        float val=tree.getFirst();int x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),k=tree.getFirstK();
        processed.add(Float.floatToIntBits(val)+","+x+","+y+","+z+","+k+"\n");tree.removeFirst();
        if(volume[k]>minimum[k]&&skel[x][y][z]==k) {
            int id=bestMembership(x,y,z),next=EMPTY;boolean direct=false;float bestPrior=-.1f;
            if(id!=k)for(int[] d:FACE) {
                int xi=x+d[0],yj=y+d[1],zl=z+d[2],neighbor=seg[xi][yj][zl];
                if(neighbor==id)direct=true;
                else if(neighbor!=k&&relation[neighbor][id]&&shape[xi][yj][zl][id]>bestPrior) {next=neighbor;bestPrior=shape[xi][yj][zl][id];}
            }
            boolean change=id!=k&&(direct||next!=EMPTY);
            if(change){if(!direct)id=next;change=topology(x,y,z,id);}
            available[x][y][z][k]=false;volume[k]-=mems[x][y][z][k];skel[x][y][z]=(byte)EMPTY;
            if(change)seg[x][y][z]=(byte)id;
            for(int[] d:FACE) {
                int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!available[xi][yj][zl][k])continue;
                int count=0;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)
                    if(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]>0&&seg[xi+i][yj+j][zl+l]==k)count++;
                float order;
                if(count==1)order=val;
                else {int best=bestMembership(xi,yj,zl);float ratio=k==best?1.0f:(mems[xi][yj][zl][best]+THIN_ZERO)/(mems[xi][yj][zl][k]+THIN_ZERO);order=val-1.0f/ratio;}
                if(volume[k]>minimum[k]&&order>-LIMIT){available[xi][yj][zl][k]=false;tree.addValue(order,xi,yj,zl,k);}
            }
        }
    }
    private void run() {float val=0.0f;while(tree.isNotEmpty()&&val>-LIMIT){val=tree.getFirst();oneEvent();}}

    private static String digest(byte[] bytes) throws Exception {
        byte[] hash=MessageDigest.getInstance("SHA-256").digest(bytes);StringBuilder out=new StringBuilder();
        for(byte b:hash)out.append(String.format("%02x",b&255));return out.toString();
    }
    private static void addFloat(MessageDigest digest,float value) {
        digest.update(ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN).putInt(Float.floatToIntBits(value)).array());
    }
    private String membershipDigest() throws Exception {
        MessageDigest d=MessageDigest.getInstance("SHA-256");for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)addFloat(d,mems[x][y][z][k]);
        return hex(d.digest());
    }
    private static String hex(byte[] hash){StringBuilder out=new StringBuilder();for(byte b:hash)out.append(String.format("%02x",b&255));return out.toString();}
    private String processedDigest() throws Exception {StringBuilder text=new StringBuilder();for(String item:processed)text.append(item);return digest(text.toString().getBytes(StandardCharsets.US_ASCII));}
    private String byteGridDigest(byte[][][] array) throws Exception {byte[] data=new byte[N*N*N];int p=0;for(byte[][] x:array)for(byte[] y:x)for(byte v:y)data[p++]=v;return digest(data);}
    private String availableDigest() throws Exception {StringBuilder out=new StringBuilder();for(boolean[][][] x:available)for(boolean[][] y:x)for(boolean[] z:y)for(boolean v:z)out.append(v?'1':'0');return digest(out.toString().getBytes(StandardCharsets.US_ASCII));}
    private String queueDigest() throws Exception {StringBuilder out=new StringBuilder();while(tree.isNotEmpty()){out.append(Float.floatToIntBits(tree.getFirst())).append(',').append(tree.getFirstX()).append(',').append(tree.getFirstY()).append(',').append(tree.getFirstZ()).append(',').append(tree.getFirstK()).append('\n');tree.removeFirst();}return digest(out.toString().getBytes(StandardCharsets.US_ASCII));}
    private static String floatBits(float[] values){StringBuilder out=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(Float.floatToIntBits(values[i]));}return out.append(']').toString();}
    private static String ints(int[] values){StringBuilder out=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]);}return out.append(']').toString();}

    public static void main(String[] args) throws Exception {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        ToadsDefaultFirstThinningReferenceHarness h=new ToadsDefaultFirstThinningReferenceHarness();
        h.initializeMembershipsAndState();h.run();
        System.out.println("{\"source\":\"VandyCRUISE MedicAlgorithmToads08 exact initial memberships -> MultipleToadSegmentation.initThinningLabels + propagateThinning\","
          +"\"preprocess_source_sha256\":\"31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7\","
          +"\"segmentation_source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\","
          +"\"algorithm_source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\","
          +"\"wrapper_source_sha256\":\"d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115\","
          +"\"lut_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\","
          +"\"heap_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\","
          +"\"shape\":["+N+","+N+","+N+","+C+"],\"crop_border\":"+BORDER+",\"limit_bits\":"+Float.floatToIntBits(LIMIT)+",\"volume_ratio_bits\":"+Float.floatToIntBits(VOLUME_RATIO)+","
          +"\"memberships_sha256\":\""+h.membershipDigest()+"\",\"initial_volume_bits\":"+floatBits(h.initialVolume)+",\"minimum_volume_bits\":"+floatBits(h.minimum)+",\"boundary_counts\":"+ints(h.boundaries)+","
          +"\"processed_count\":"+h.processed.size()+",\"processed_sha256\":\""+h.processedDigest()+"\",\"segmentation_sha256\":\""+h.byteGridDigest(h.seg)+"\",\"skeleton_sha256\":\""+h.byteGridDigest(h.skel)+"\",\"available_sha256\":\""+h.availableDigest()+"\",\"final_volume_bits\":"+floatBits(h.volume)+",\"queue_after_count\":0,\"queue_after_sha256\":\""+h.queueDigest()+"\"}");
    }
}
