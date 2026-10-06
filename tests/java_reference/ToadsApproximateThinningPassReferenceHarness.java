import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;

/** Small source-order executable for the sparse approximate queue pass. */
public final class ToadsApproximateThinningPassReferenceHarness {
    static final int N=9,C=3,A=2;
    static final class Heap {
        float[] v=new float[N*N*N*4+1]; int[] x=new int[v.length],y=new int[v.length],z=new int[v.length],k=new int[v.length]; int size=0;
        void add(float val,int a,int b,int c,int d){int h=++size;while(val>v[h/2]){v[h]=v[h/2];x[h]=x[h/2];y[h]=y[h/2];z[h]=z[h/2];k[h]=k[h/2];h/=2;}v[h]=val;x[h]=a;y[h]=b;z[h]=c;k[h]=d;}
        void remove(){int h=1,child;float val=v[size];int a=x[size],b=y[size],c=z[size],d=k[size];size--;for(;h*2<=size;h=child){child=h*2;if(child!=size&&v[child+1]>v[child])child++;if(v[child]>val){v[h]=v[child];x[h]=x[child];y[h]=y[child];z[h]=z[child];k[h]=k[child];}else break;}v[h]=val;x[h]=a;y[h]=b;z[h]=c;k[h]=d;}
    }
    static int cls(int x,int y,int z){int fold=((x*3+y*2)%3)-1,dx=x-4,dy=y-4,dz=z-4,r=2+fold;if(dx*dx+dy*dy+dz*dz<r*r)return 1;int fz=4+((x+2*y)%3)-1;return Math.abs(z-fz)<=1&&Math.abs(y-4)<=3&&Math.abs(x-4)<=3?2:0;}
    static int top(int x,int y,int z,int a){return a==0?cls(x,y,z):(cls(x,y,z)+1)%C;}
    static float mem(int x,int y,int z,int a){return (((x*17+y*11+z*7+a*13)%39)+1)/100.0f;}
    static int count(int[][][] s,int x,int y,int z,int k){int n=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(s[x+i][y+j][z+l]==k)n++;return n;}
    static String sha(byte[] value)throws Exception{byte[] d=MessageDigest.getInstance("SHA-256").digest(value);StringBuilder s=new StringBuilder();for(byte b:d)s.append(String.format("%02x",b));return s.toString();}
    public static void main(String[] args)throws Exception{
        int[][][] seg=new int[N][N][N],skel=new int[N][N][N];boolean[][][][] av=new boolean[N][N][N][C];float[] vol=new float[C];Heap q=new Heap();
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++){seg[x][y][z]=cls(x,y,z);skel[x][y][z]=-1;}
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){int cur=seg[x][y][z];skel[x][y][z]=cur;for(int a=0;a<A;a++)if(cur==top(x,y,z,a))vol[cur]+=mem(x,y,z,a);for(int k=0;k<C;k++)av[x][y][z][k]=k==cur;}
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++){int cur=seg[x][y][z];boolean boundary=seg[x-1][y][z]!=cur||seg[x+1][y][z]!=cur||seg[x][y-1][z]!=cur||seg[x][y+1][z]!=cur||seg[x][y][z-1]!=cur||seg[x][y][z+1]!=cur;if(!boundary)continue;av[x][y][z][cur]=false;q.add(-1.0f,x,y,z,cur);}
        MessageDigest events=MessageDigest.getInstance("SHA-256");int processed=0;float val=0;final float limit=1000.0f;
        while(q.size>0&&val>-limit){int x=q.x[1],y=q.y[1],z=q.z[1],k=q.k[1];val=q.v[1];q.remove();ByteBuffer b=ByteBuffer.allocate(20).order(ByteOrder.BIG_ENDIAN).putInt(Float.floatToIntBits(val)).putInt(x).putInt(y).putInt(z).putInt(k);events.update(b.array());processed++;
            if(vol[k]>0.0f&&skel[x][y][z]==k){av[x][y][z][k]=false;for(int a=0;a<A;a++)if(top(x,y,z,a)==k)vol[k]-=mem(x,y,z,a);skel[x][y][z]=-1;
                for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l==1){int xi=x+i,yj=y+j,zl=z+l;if(!av[xi][yj][zl][k])continue;int nn=count(seg,xi,yj,zl,k);float order;if(nn==1)order=val;else{int t=top(xi,yj,zl,0);float ratio;if(k==t)ratio=1.0f;else{float cur=0.0f;for(int a=1;a<A;a++)if(top(xi,yj,zl,a)==k)cur=mem(xi,yj,zl,a);ratio=(mem(xi,yj,zl,0)+1e-10f)/(cur+1e-10f);}order=val-1.0f/ratio;}if(vol[k]>0.0f&&order>-limit){av[xi][yj][zl][k]=false;q.add(order,xi,yj,zl,k);}}
            }
        }
        System.out.printf("{\"processed\":%d,\"remaining\":%d,\"event_sha256\":\"%s\",\"volume_bits\":[%d,%d,%d],\"skeleton_remaining\":%d}%n",processed,q.size,hex(events.digest()),Float.floatToIntBits(vol[0]),Float.floatToIntBits(vol[1]),Float.floatToIntBits(vol[2]),countSkeleton(skel));
    }
    static int countSkeleton(int[][][] s){int n=0;for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)if(s[x][y][z]>=0)n++;return n;}
    static String hex(byte[] d){StringBuilder s=new StringBuilder();for(byte b:d)s.append(String.format("%02x",b));return s.toString();}
}
