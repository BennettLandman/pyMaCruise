/** Source-ordered six-face rescheduling excerpt from propagateThinning. */
public final class ToadsApproximateNeighborScheduleReferenceHarness {
    static final int N=9, C=3, A=2;
    static int cls(int x,int y,int z) {
        int folded=((x*3+y*2)%3)-1, dx=x-4,dy=y-4,dz=z-4;
        int radius=2+folded;
        if(dx*dx+dy*dy+dz*dz<radius*radius) return 1;
        int foldZ=4+((x+2*y)%3)-1;
        return Math.abs(z-foldZ)<=1 && Math.abs(y-4)<=3 && Math.abs(x-4)<=3 ? 2 : 0;
    }
    static int top(int x,int y,int z,int slot) { return slot==0 ? cls(x,y,z) : (cls(x,y,z)+1)%C; }
    static float mem(int x,int y,int z,int slot) { return (((x*17+y*11+z*7+slot*13)%39)+1)/100.0f; }
    static int count(int[][][] seg,int x,int y,int z,int k) {
        int n=0; for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++) if(seg[x+i][y+j][z+l]==k)n++;
        return n;
    }
    public static void main(String[] args) {
        int[][][] seg=new int[N][N][N]; boolean[][][][] available=new boolean[N][N][N][C];
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)seg[x][y][z]=cls(x,y,z);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<C;k++)available[x][y][z][k]=seg[x][y][z]==k;
        // The pinned first one-step event is (1,1,1), class 0, priority -1.
        int x=1,y=1,z=1,k=0; float val=-1.0f; StringBuilder out=new StringBuilder("{\"inserted\":["); boolean first=true;
        for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l==1){
            int xi=x+i,yj=y+j,zl=z+l; if(!available[xi][yj][zl][k])continue;
            int neighborCount=count(seg,xi,yj,zl,k); float order;
            if(neighborCount==1)order=val;
            else {
                int top=top(xi,yj,zl,0); float ratio;
                if(k==top)ratio=1.0f;
                else {float cur=0.0f;for(int a=1;a<A;a++)if(top(xi,yj,zl,a)==k)cur=mem(xi,yj,zl,a); ratio=(mem(xi,yj,zl,0)+1e-10f)/(cur+1e-10f);}
                order=val-1.0f/ratio;
            }
            if(1000.0f>0.0f && order>-1000.0f){
                available[xi][yj][zl][k]=false;
                if(!first)out.append(','); first=false;
                out.append(String.format("[%d,%d,%d,%d]",Float.floatToIntBits(order),xi,yj,zl));
            }
        }
        out.append("]}"); System.out.println(out);
    }
}
