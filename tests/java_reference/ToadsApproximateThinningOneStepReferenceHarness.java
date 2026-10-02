/** One selected approximate propagateThinning event after initThinningLabels. */
public final class ToadsApproximateThinningOneStepReferenceHarness {
    static final int N=9, C=3, A=2;
    static final class Heap {
        float[] value = new float[N*N*N+1];
        int[] x = new int[N*N*N+1], y = new int[N*N*N+1], z = new int[N*N*N+1], label = new int[N*N*N+1];
        int size=0;
        void add(float v,int a,int b,int c,int k) {
            int hole=++size;
            for (; v>value[hole/2]; hole/=2) {
                value[hole]=value[hole/2]; x[hole]=x[hole/2]; y[hole]=y[hole/2]; z[hole]=z[hole/2]; label[hole]=label[hole/2];
            }
            value[hole]=v; x[hole]=a; y[hole]=b; z[hole]=c; label[hole]=k;
        }
        void remove() {
            int hole=1; float v=value[size]; int a=x[size],b=y[size],c=z[size],k=label[size]; size--;
            int child;
            for (; hole*2<=size; hole=child) {
                child=hole*2;
                if (child!=size && value[child+1]>value[child]) child++;
                if (value[child]>v) { value[hole]=value[child]; x[hole]=x[child]; y[hole]=y[child]; z[hole]=z[child]; label[hole]=label[child]; }
                else break;
            }
            value[hole]=v; x[hole]=a; y[hole]=b; z[hole]=c; label[hole]=k;
        }
    }
    static int cls(int x,int y,int z) {
        int folded=((x*3+y*2)%3)-1;
        int dx=x-4,dy=y-4,dz=z-4;
        int radius=2+folded;
        int radius2=dx*dx+dy*dy+dz*dz;
        if (radius2 < radius*radius) return 1;
        int foldZ=4+((x+2*y)%3)-1;
        if (Math.abs(z-foldZ)<=1 && Math.abs(y-4)<=3 && Math.abs(x-4)<=3) return 2;
        return 0;
    }
    static int top(int x,int y,int z,int slot) { return slot==0 ? cls(x,y,z) : (cls(x,y,z)+1)%C; }
    static float mem(int x,int y,int z,int slot) {
        int numerator=(x*17+y*11+z*7+slot*13)%39+1;
        return numerator/100.0f;
    }
    static int bits(float value) { return Float.floatToIntBits(value); }
    public static void main(String[] args) {
        int[][][] seg=new int[N][N][N], skel=new int[N][N][N];
        boolean[][][][] available=new boolean[N][N][N][C];
        float[] volume=new float[C];
        for(int x=0;x<N;x++) for(int y=0;y<N;y++) for(int z=0;z<N;z++) { seg[x][y][z]=cls(x,y,z); skel[x][y][z]=-1; }
        for(int x=1;x<N-1;x++) for(int y=1;y<N-1;y++) for(int z=1;z<N-1;z++) {
            int current=seg[x][y][z]; skel[x][y][z]=current;
            for(int a=0;a<A;a++) if(current==top(x,y,z,a)) volume[current]+=mem(x,y,z,a);
            for(int k=0;k<C;k++) available[x][y][z][k]=(current==k);
        }
        Heap queue=new Heap();
        for(int x=1;x<N-1;x++) for(int y=1;y<N-1;y++) for(int z=1;z<N-1;z++) {
            int current=seg[x][y][z]; boolean boundary=false; int neighborCount=0,ni=0,nj=0,nl=0;
            for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++) {
                int d2=i*i+j*j+l*l;
                if(d2==1 && seg[x+i][y+j][z+l]!=current) boundary=true;
                else if(d2>0) { neighborCount++; ni=x+i; nj=y+j; nl=z+l; }
            }
            if(!boundary) continue;
            available[x][y][z][current]=false;
            boolean endpoint=false;
            if(neighborCount==1) {
                neighborCount=0;
                for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++)
                    if(i*i+j*j+l*l>0 && seg[ni+i][nj+j][nl+l]==current) neighborCount++;
                endpoint=neighborCount<=2;
            }
            float priority=endpoint ? -0.0f : -1.0f;
            queue.add(priority,x,y,z,current);
        }
        int x=queue.x[1],y=queue.y[1],z=queue.z[1],k=queue.label[1]; float priority=queue.value[1];
        queue.remove();
        // best[...,0] is constructed equal to segmentation, so source id==k:
        // no label change, but the eligible voxel is removed from the skeleton.
        available[x][y][z][k]=false;
        float decrement=0.0f;
        for(int a=0;a<A;a++) if(top(x,y,z,a)==k) decrement+=mem(x,y,z,a);
        volume[k]-=decrement;
        skel[x][y][z]=-1;
        System.out.printf("{\"shape\":[9,9,9],\"popped\":[%d,%d,%d,%d,%d],\"volume_bits\":[%d,%d,%d],\"center_label\":%d,\"center_skeleton\":%d,\"center_available\":[%s,%s,%s],\"queue_remaining\":%d}%n",
            bits(priority),x,y,z,k,bits(volume[0]),bits(volume[1]),bits(volume[2]),seg[x][y][z],skel[x][y][z],
            available[x][y][z][0],available[x][y][z][1],available[x][y][z][2],queue.size);
    }
}
