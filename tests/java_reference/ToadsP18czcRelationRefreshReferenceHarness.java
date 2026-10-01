import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;

/** Selected VandyCRUISE relation refresh on caller-supplied post-growth arrays. */
public final class ToadsP18czcRelationRefreshReferenceHarness {
    private static final String SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static byte[] readBytes(Path path, int count) throws Exception {
        byte[] values = Files.readAllBytes(path);
        if (values.length != count) throw new IllegalArgumentException("unexpected byte payload: " + path);
        return values;
    }
    private static float[] readFloats(Path path, int count) throws Exception {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length != count * 4) throw new IllegalArgumentException("unexpected float payload: " + path);
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[count];
        for (int i = 0; i < count; i++) values[i] = buffer.getFloat();
        return values;
    }
    private static byte[][][] labels(byte[] flat, int nx, int ny, int nz) {
        byte[][][] values = new byte[nx][ny][nz]; int p = 0;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) values[x][y][z]=flat[p++];
        return values;
    }
    private static float[][][][] memberships(float[] flat, int nx, int ny, int nz, int classes) {
        float[][][][] values = new float[nx][ny][nz][classes]; int p = 0;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
            for (int k=0;k<classes;k++) values[x][y][z][k]=flat[p++];
        return values;
    }
    private static byte[][][][] boundary(byte[] flat, int nx, int ny, int nz, int classes) {
        byte[][][][] values = new byte[nx][ny][nz][classes]; int p = 0;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
            for (int k=0;k<classes;k++) values[x][y][z][k]=flat[p++];
        return values;
    }
    private static void append(StringBuilder out, int[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if(i>0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    private static int[] flatten(short[][][] values) {
        int[] out=new int[values.length*values[0].length*values[0][0].length]; int p=0;
        for(short[][] x:values) for(short[] y:x) for(short z:y) out[p++]=z;
        return out;
    }
    private static int[] flatten(byte[][][][] values) {
        int[] out=new int[values.length*values[0].length*values[0][0].length*values[0][0][0].length]; int p=0;
        for(byte[][][] x:values) for(byte[][] y:x) for(byte[] z:y) for(byte k:z) out[p++]=k;
        return out;
    }
    private static int[] flatten(byte[][][][] values, boolean closest) {
        int[] out=new int[values.length*values[0].length*values[0][0].length*(closest?4:values[0][0][0].length)]; int p=0;
        for(byte[][][] x:values) for(byte[][] y:x) for(byte[] z:y) for(byte k:z) out[p++]=k;
        return out;
    }
    private static int[] flatten(boolean[][] values) {
        int[] out=new int[values.length*values[0].length]; int p=0;
        for(boolean[] x:values) for(boolean y:x) out[p++]=y?1:0;
        return out;
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=8) throw new IllegalArgumentException("seg.bin mem.f32le boundary.bin nx ny nz classes limit required");
        int nx=Integer.parseInt(args[3]), ny=Integer.parseInt(args[4]), nz=Integer.parseInt(args[5]);
        int classes=Integer.parseInt(args[6]), limit=Integer.parseInt(args[7]);
        byte[][][] segmentation=labels(readBytes(Path.of(args[0]),nx*ny*nz),nx,ny,nz);
        float[][][][] mem=memberships(readFloats(Path.of(args[1]),nx*ny*nz*classes),nx,ny,nz,classes);
        byte[][][][] initial=boundary(readBytes(Path.of(args[2]),nx*ny*nz*classes),nx,ny,nz,classes);
        ToadsRelationListReferenceHarness.Result list=ToadsRelationListReferenceHarness.build(segmentation,classes);
        ToadsRelationDiffusionReferenceHarness.Outcome relation=ToadsRelationDiffusionReferenceHarness.computeRelations(segmentation,classes,mem,limit);
        ToadsBoundaryDistancesReferenceHarness.Outcome distances=ToadsBoundaryDistancesReferenceHarness.compute(segmentation,initial,classes,limit);
        byte[][][][] closest=ToadsClosestLabelsReferenceHarness.closestLabels(segmentation,limit);
        StringBuilder table=new StringBuilder("[");
        for(int i=0;i<list.relations.length;i++){
            if(i>0)table.append(','); table.append('[');
            for(int k=0;k<classes;k++){if(k>0)table.append(',');table.append(list.relations[i][k]?1:0);}
            table.append(']');
        }
        System.out.println("{\"source_sha256\":\""+SOURCE_SHA+"\",\"relation_diffusion_iterations\":"+relation.iterations
            +",\"boundary_diffusion_passes\":"+distances.passes+",\"relation_map\":"+array(flatten(relation.map))
            +",\"relation_membership\":"+table+"],\"boundary_distances\":"+array(flatten(distances.distances))
            +",\"closest_labels\":"+array(flatten(closest,true))+"}");
    }
    private static String array(int[] values){StringBuilder out=new StringBuilder();append(out,values);return out.toString();}
}
