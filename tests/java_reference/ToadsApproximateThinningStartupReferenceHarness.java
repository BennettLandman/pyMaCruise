import java.util.ArrayList;
import java.util.List;

/** Selected useApprox arm of VandyCRUISE MultipleToadSegmentation.initThinningLabels. */
public final class ToadsApproximateThinningStartupReferenceHarness {
    static final int N = 5, C = 3, A = 2;
    static final float ZERO = 1.0e-30f;

    static final class Event {
        float value; int x, y, z, label;
        Event(float value, int x, int y, int z, int label) {
            this.value = value; this.x = x; this.y = y; this.z = z; this.label = label;
        }
    }

    static int label(int x, int y, int z) { return (x + 2 * y + z) % C; }
    static int best(int x, int y, int z, int a) { return (x + y + z + a) % C; }
    static float membership(int x, int y, int z, int a) {
        return ((x * 17 + y * 11 + z * 7 + a * 13) % 79 + 1) / 100.0f;
    }
    static int bits(float value) { return Float.floatToIntBits(value); }

    static void append(StringBuilder out, String value) {
        if (out.length() > 0) out.append(',');
        out.append(value);
    }

    public static void main(String[] args) {
        int[][][] segmentation = new int[N][N][N];
        int[][][] skeleton = new int[N][N][N];
        boolean[][][][] available = new boolean[N][N][N][C];
        float[] volumes = new float[C];
        int[] boundaries = new int[C];
        for (int x=0; x<N; x++) for (int y=0; y<N; y++) for (int z=0; z<N; z++) {
            segmentation[x][y][z] = label(x,y,z);
            skeleton[x][y][z] = -9;
            for (int c=0; c<C; c++) available[x][y][z][c] = (x==0 || y==0 || z==0 || x==N-1 || y==N-1 || z==N-1) && c==1;
        }
        List<Event> events = new ArrayList<>();
        for (int x=1; x<N-1; x++) for (int y=1; y<N-1; y++) for (int z=1; z<N-1; z++) {
            int current = segmentation[x][y][z];
            skeleton[x][y][z] = current;
            for (int a=0; a<A; a++) if (current == best(x,y,z,a)) volumes[current] += membership(x,y,z,a);
            for (int c=0; c<C; c++) available[x][y][z][c] = c==current;
        }
        for (int x=1; x<N-1; x++) for (int y=1; y<N-1; y++) for (int z=1; z<N-1; z++) {
            int current = segmentation[x][y][z];
            boolean boundary = false; int neighborCount = 0, ni=0, nj=0, nl=0;
            for (int i=-1; i<=1; i++) for (int j=-1; j<=1; j++) for (int k=-1; k<=1; k++) {
                int d2=i*i+j*j+k*k, neighbor=segmentation[x+i][y+j][z+k];
                if (d2==1 && neighbor!=current) boundary=true;
                else if (d2>0) { neighborCount++; ni=x+i; nj=y+j; nl=z+k; }
            }
            if (!boundary) continue;
            available[x][y][z][current]=false;
            boolean endpoint=false;
            if (neighborCount==1) {
                neighborCount=0;
                for (int i=-1; i<=1; i++) for (int j=-1; j<=1; j++) for (int k=-1; k<=1; k++)
                    if (i*i+j*j+k*k>0 && segmentation[ni+i][nj+j][nl+k]==current) neighborCount++;
                endpoint=neighborCount<=2;
            }
            float priority;
            if (endpoint) priority=-0.0f;
            else {
                int top=best(x,y,z,0); float ratio;
                if (current==top) ratio=1.0f;
                else {
                    float cur=0.0f;
                    for (int a=1; a<A; a++) if (best(x,y,z,a)==current) cur=membership(x,y,z,a);
                    ratio=(membership(x,y,z,0)+ZERO)/(cur+ZERO);
                }
                priority=-1.0f/ratio;
            }
            boundaries[current]++;
            events.add(new Event(priority,x,y,z,current));
        }
        StringBuilder out=new StringBuilder("{\"shape\":[5,5,5],\"classes\":3,\"slots\":2,\"volumes_bits\":[");
        for (int c=0;c<C;c++) { if(c>0) out.append(','); out.append(bits(volumes[c])); }
        out.append("],\"boundary_counts\":[");
        for (int c=0;c<C;c++) { if(c>0) out.append(','); out.append(boundaries[c]); }
        out.append("],\"skeleton\":\"");
        for (int x=0;x<N;x++) for (int y=0;y<N;y++) for (int z=0;z<N;z++) out.append(skeleton[x][y][z]);
        out.append("\",\"available\":\"");
        for (int x=0;x<N;x++) for (int y=0;y<N;y++) for (int z=0;z<N;z++) for (int c=0;c<C;c++) out.append(available[x][y][z][c]?'1':'0');
        out.append("\",\"priority_events\":[");
        for (int index=0; index<events.size(); index++) {
            Event e=events.get(index); if(index>0) out.append(',');
            out.append('[').append(bits(e.value)).append(',').append(e.x).append(',').append(e.y).append(',').append(e.z).append(',').append(e.label).append(']');
        }
        out.append("]}");
        System.out.println(out);
    }
}
