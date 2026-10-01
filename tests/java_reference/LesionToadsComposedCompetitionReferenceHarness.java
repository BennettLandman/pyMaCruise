import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/** Selected competition excerpts applied to the preceding folded synthetic continuation. */
public final class LesionToadsComposedCompetitionReferenceHarness {
    private static final int CLASSES = 13;
    private static final String[] TYPES = {
        "object", "object", "object", "object", "object", "object", "object",
        "object", "object", "object", "lesion", "object", "object"
    };
    private static float[] initialVolume;
    private static CriticalPointLUT lut;

    private static final class State {
        int nx, ny, nz, channels;
        short[][][] segmentation;
        float[][][][] mems;
        float[][][] lesion, blackhole, ordering;
    }

    private static State readNifti(String path) throws Exception {
        byte[] bytes;
        try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(path));
                ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            in.transferTo(out);
            bytes = out.toByteArray();
        }
        ByteBuffer b = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        if (b.getInt(0) != 348 || b.getShort(70) != 16 || b.getShort(72) != 32)
            throw new IllegalArgumentException("expected little-endian NIfTI-1 float32");
        int dims = b.getShort(40);
        State s = new State();
        s.nx = b.getShort(42); s.ny = b.getShort(44); s.nz = b.getShort(46);
        s.channels = dims >= 4 ? b.getShort(48) : 1;
        if (s.channels != CLASSES + 4) throw new IllegalArgumentException("unexpected frames");
        int offset = (int)b.getFloat(108);
        float[] raw = new float[s.nx * s.ny * s.nz * s.channels];
        for (int i=0;i<raw.length;i++) raw[i] = b.getFloat(offset + 4*i);
        s.segmentation = new short[s.nx][s.ny][s.nz];
        s.mems = new float[s.nx][s.ny][s.nz][CLASSES];
        s.lesion = new float[s.nx][s.ny][s.nz];
        s.blackhole = new float[s.nx][s.ny][s.nz];
        s.ordering = new float[s.nx][s.ny][s.nz];
        for (int x=0;x<s.nx;x++) for (int y=0;y<s.ny;y++) for (int z=0;z<s.nz;z++) {
            s.segmentation[x][y][z] = (short)raw[index(s,x,y,z,0)];
            for (int k=0;k<CLASSES;k++) s.mems[x][y][z][k] = raw[index(s,x,y,z,1+k)];
            s.lesion[x][y][z] = raw[index(s,x,y,z,1+CLASSES)];
            s.blackhole[x][y][z] = raw[index(s,x,y,z,2+CLASSES)];
            s.ordering[x][y][z] = raw[index(s,x,y,z,3+CLASSES)];
        }
        return s;
    }

    private static int index(State s, int x, int y, int z, int c) {
        return x + s.nx*(y + s.ny*(z + s.nz*c));
    }

    private static State copy(State source) {
        State target = new State();
        target.nx=source.nx; target.ny=source.ny; target.nz=source.nz; target.channels=source.channels;
        target.segmentation=new short[target.nx][target.ny][target.nz];
        target.mems=new float[target.nx][target.ny][target.nz][CLASSES];
        target.lesion=new float[target.nx][target.ny][target.nz];
        target.blackhole=new float[target.nx][target.ny][target.nz];
        target.ordering=new float[target.nx][target.ny][target.nz];
        for(int x=0;x<target.nx;x++) for(int y=0;y<target.ny;y++) for(int z=0;z<target.nz;z++) {
            target.segmentation[x][y][z]=source.segmentation[x][y][z];
            System.arraycopy(source.mems[x][y][z],0,target.mems[x][y][z],0,CLASSES);
            target.lesion[x][y][z]=source.lesion[x][y][z];
            target.blackhole[x][y][z]=source.blackhole[x][y][z];
            target.ordering[x][y][z]=source.ordering[x][y][z];
        }
        return target;
    }

    private static boolean homeomorphic(short[][][] seg, int x, int y, int z, short k) {
        if (seg[x][y][z] == k) return true;
        boolean[][][] obj = new boolean[3][3][3];
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
            obj[1+i][1+j][1+l] = seg[x+i][y+j][z+l] == k;
        obj[1][1][1] = true;
        if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
            obj[1+i][1+j][1+l] = seg[x+i][y+j][z+l] == seg[x][y][z];
        obj[1][1][1] = false;
        if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        int count=0;
        short[] labels = new short[26];
        for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
            if ((i*i+j*j+l*l>0) && seg[x+i][y+j][z+l] != k
                    && seg[x+i][y+j][z+l] != seg[x][y][z]) {
                boolean found=false;
                for (int n=0;n<count;n++) if (seg[x+i][y+j][z+l] == labels[n]) { found=true; break; }
                if (!found) labels[count++] = seg[x+i][y+j][z+l];
            }
        }
        for (int n=0;n<count;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
                obj[1+i][1+j][1+l] = seg[x+i][y+j][z+l] == seg[x][y][z]
                        || seg[x+i][y+j][z+l] == labels[n];
            obj[1][1][1]=false;
            if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
                obj[1+i][1+j][1+l] = seg[x+i][y+j][z+l] == k
                        || seg[x+i][y+j][z+l] == labels[n];
            obj[1][1][1]=true;
            if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) for (int m=n+1;m<count;m++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                short label=seg[x+i][y+j][z+l];
                obj[1+i][1+j][1+l] = label == seg[x][y][z] || label == labels[n] || label == labels[m];
            }
            obj[1][1][1]=false;
            if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        }
        for (int n=0;n<count;n++) for (int m=n+1;m<count;m++) {
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
                short label=seg[x+i][y+j][z+l];
                obj[1+i][1+j][1+l] = label == k || label == labels[n] || label == labels[m];
            }
            obj[1][1][1]=true;
            if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
        }
        return true;
    }

    private static float initialOrder(State s, int x, int y, int z, int old, int target, int lesionClass) {
        if (target == lesionClass)
            return (s.mems[x][y][z][target]+s.lesion[x][y][z]+s.blackhole[x][y][z])
                    /(s.mems[x][y][z][old]+1e-10f);
        if (old == lesionClass)
            return s.mems[x][y][z][target]
                    /(s.mems[x][y][z][old]+s.lesion[x][y][z]+s.blackhole[x][y][z]+1e-10f);
        return s.mems[x][y][z][target]/(s.mems[x][y][z][old]+1e-10f);
    }

    private static Float neighborOrder(State s, int x, int y, int z, int old, int target, int lesionClass) {
        if (old == lesionClass) {
            if (TYPES[old].equals("lesion"))
                return s.mems[x][y][z][target]/(s.mems[x][y][z][old]+s.lesion[x][y][z]+1e-10f);
            if (TYPES[old].equals("blackhole"))
                return s.mems[x][y][z][target]/(s.mems[x][y][z][old]+s.lesion[x][y][z]+s.blackhole[x][y][z]+1e-10f);
        } else if (target == lesionClass) {
            if (TYPES[old].equals("lesion"))
                return (s.mems[x][y][z][target]+s.lesion[x][y][z])/(s.mems[x][y][z][old]+1e-10f);
            if (TYPES[old].equals("blackhole"))
                return (s.mems[x][y][z][target]+s.lesion[x][y][z]+s.blackhole[x][y][z])/(s.mems[x][y][z][old]+1e-10f);
        } else return s.mems[x][y][z][target]/(s.mems[x][y][z][old]+1e-10f);
        return null;
    }

    private static int bits(float f) { return Float.floatToRawIntBits(f); }
    private static String hex(byte[] bytes) {
        StringBuilder s=new StringBuilder(); for(byte value:bytes) s.append(String.format("%02x",value)); return s.toString();
    }
    private static String digest(short[][][] a) throws Exception {
        ByteBuffer b=ByteBuffer.allocate(a.length*a[0].length*a[0][0].length*2).order(ByteOrder.LITTLE_ENDIAN);
        for(short[][] x:a) for(short[] y:x) for(short v:y) b.putShort(v);
        return hex(MessageDigest.getInstance("SHA-256").digest(b.array()));
    }
    private static String digest(float[][][] a) throws Exception {
        ByteBuffer b=ByteBuffer.allocate(a.length*a[0].length*a[0][0].length*4).order(ByteOrder.LITTLE_ENDIAN);
        for(float[][] x:a) for(float[] y:x) for(float v:y) b.putInt(bits(v));
        return hex(MessageDigest.getInstance("SHA-256").digest(b.array()));
    }
    private static String digest(boolean[][][][] a) throws Exception {
        ByteArrayOutputStream b=new ByteArrayOutputStream(); for(boolean[][][] x:a) for(boolean[][] y:x) for(boolean[] z:y) for(boolean v:z) b.write(v?1:0);
        return hex(MessageDigest.getInstance("SHA-256").digest(b.toByteArray()));
    }
    private static String digest(float[] a) throws Exception {
        ByteBuffer b=ByteBuffer.allocate(a.length*4).order(ByteOrder.LITTLE_ENDIAN); for(float v:a) b.putInt(bits(v));
        return hex(MessageDigest.getInstance("SHA-256").digest(b.array()));
    }
    private static String item(float value, int x, int y, int z, int payload) {
        return "["+bits(value)+","+x+","+y+","+z+","+payload+"]";
    }

    private static String run(State s, int lesionClass) throws Exception {
        short[][][] seg=s.segmentation;
        float[] volume=initialVolume.clone();
        boolean[][][][] available=new boolean[s.nx][s.ny][s.nz][CLASSES];
        short[][][] skeleton=new short[s.nx][s.ny][s.nz];
        BinaryHeap4D tree=new BinaryHeap4D(s.nx*s.ny*s.nz*6, BinaryHeap4D.MAXTREE);
        tree.reset(); tree.setMaxTree();
        List<String> initial=new ArrayList<>();
        for(int x=0;x<s.nx;x++) for(int y=0;y<s.ny;y++) for(int z=0;z<s.nz;z++) {
            for(int c=0;c<CLASSES;c++) available[x][y][z][c]=true;
            skeleton[x][y][z]=seg[x][y][z];
        }
        for(int x=1;x<s.nx-1;x++) for(int y=1;y<s.ny-1;y++) for(int z=1;z<s.nz-1;z++)
            for(byte k=0;k<CLASSES;k++) if(seg[x][y][z]==k)
                for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int l=-1;l<=1;l++) if(i*i+j*j+l*l==1) {
                    int xi=x+i,yj=y+j,zl=z+l; int target=seg[xi][yj][zl];
                    if(target!=k) {
                        float order=initialOrder(s,xi,yj,zl,k,target,lesionClass);
                        if(order>0) { int packed=k+CLASSES*target; tree.addValue(order,x,y,z,packed); initial.add(item(order,x,y,z,packed)); }
                    }
                }
        StringBuilder events=new StringBuilder("[");
        int eventCount=0, acceptedCount=0, rejectedCount=0, staleCount=0;
        float order=0.0f;
        while(tree.isNotEmpty()) {
            short x=tree.getFirstX(),y=tree.getFirstY(),z=tree.getFirstZ(),packed=tree.getFirstK();
            float value=tree.getFirst(); tree.removeFirst();
            int target=packed/CLASSES, old=packed-target*CLASSES;
            boolean processed=seg[x][y][z]==old; Boolean accepted=null;
            List<String> inserted=new ArrayList<>();
            if(processed) {
                accepted=homeomorphic(seg,x,y,z,(short)target);
                if(accepted) {
                    seg[x][y][z]=(short)target; skeleton[x][y][z]=(short)target; s.ordering[x][y][z]=value;
                    if(target==lesionClass) {
                        volume[target]+=s.mems[x][y][z][target]+s.lesion[x][y][z];
                        if(TYPES[target].equals("blackhole")) volume[target]+=s.blackhole[x][y][z];
                    } else volume[target]+=s.mems[x][y][z][target];
                    if(old==lesionClass) {
                        volume[old]-=s.mems[x][y][z][old]-s.lesion[x][y][z];
                        if(TYPES[old].equals("blackhole")) volume[old]-=s.blackhole[x][y][z];
                    } else volume[old]-=s.mems[x][y][z][old];
                    int[][] faces={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
                    for(int[] d:faces) {
                        int xi=x+d[0],yj=y+d[1],zl=z+d[2]; if(seg[xi][yj][zl]==old) {
                            Float computed=neighborOrder(s,xi,yj,zl,old,target,lesionClass);
                            if(computed!=null) order=computed;
                            if(order>0) { tree.addValue(order,xi,yj,zl,old); inserted.add(item(order,xi,yj,zl,old)); }
                        }
                    }
                    acceptedCount++;
                } else rejectedCount++;
            } else staleCount++;
            if(eventCount++>0) events.append(',');
            events.append("{\"pop\":").append(item(value,x,y,z,packed)).append(",\"processed\":").append(processed)
                    .append(",\"candidate\":").append(target).append(",\"accepted\":").append(accepted==null?"null":accepted)
                    .append(",\"inserted\":[").append(String.join(",",inserted)).append("]}");
        }
        events.append(']');
        int[] volumeBits=new int[volume.length]; for(int i=0;i<volume.length;i++) volumeBits[i]=bits(volume[i]);
        return "{\"lesion_class\":"+lesionClass+",\"initial_queue\":["+String.join(",",initial)+"],\"events\":"+events
                +",\"accepted\":"+acceptedCount+",\"rejected\":"+rejectedCount+",\"stale\":"+staleCount
                +",\"final_segmentation_sha256\":\""+digest(seg)+"\",\"final_skeleton_sha256\":\""+digest(skeleton)
                +"\",\"availability_sha256\":\""+digest(available)+"\",\"ordering_sha256\":\""+digest(s.ordering)
                +"\",\"volume_bits\":"+java.util.Arrays.toString(volumeBits)+",\"remaining_queue\":0}";
    }

    public static void main(String[] args) throws Exception {
        PrintStream output=System.out; System.setOut(new PrintStream(new ByteArrayOutputStream()));
        lut=new CriticalPointLUT("/edu/jhmi/rad/medic/structures/critical186LUT.raw.gz",200);
        boolean loaded=lut.loadCompressedPattern(); System.setOut(output);
        if(!loaded) throw new IllegalStateException("unable to load selected 18/6 LUT");
        if (args.length != 2) throw new IllegalArgumentException("expected NIfTI path and comma-separated initial volume bits");
        String[] volumeBits=args[1].split(",");
        if (volumeBits.length != CLASSES) throw new IllegalArgumentException("expected one volume per class");
        initialVolume=new float[CLASSES];
        for(int i=0;i<CLASSES;i++) initialVolume[i]=Float.intBitsToFloat((int)Long.parseLong(volumeBits[i],16));
        State source=readNifti(args[0]);
        System.out.println(run(copy(source),10));
    }
}
