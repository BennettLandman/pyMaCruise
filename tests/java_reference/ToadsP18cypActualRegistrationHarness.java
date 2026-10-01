import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.lang.reflect.Field;
import edu.jhmi.rad.medic.methods.MultipleToadDeformableAtlas;
import edu.jhmi.rad.medic.methods.ParametricTransform;

/** Calls the supplied VandyCRUISE MultipleToadDeformableAtlas registration methods headlessly. */
public final class ToadsP18cypActualRegistrationHarness {
    private static float[] read(Path path, int count) throws Exception {
        byte[] bytes = Files.readAllBytes(path);
        if (bytes.length != count * 4) throw new IllegalArgumentException("unexpected float payload: " + path);
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[count];
        for (int i=0; i<count; i++) values[i] = buffer.getFloat();
        return values;
    }
    private static void set(Object target, String name, Object value) throws Exception {
        Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        field.set(target, value);
    }
    private static Object get(Object target, String name) throws Exception {
        Field field = target.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.get(target);
    }
    private static String bits(float value) { return Integer.toString(Float.floatToIntBits(value)); }
    private static String vector(float[] values) {
        StringBuilder out=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(bits(values[i]));}
        return out.append(']').toString();
    }
    private static String bounds(MultipleToadDeformableAtlas atlas, int classes) throws Exception {
        StringBuilder out=new StringBuilder("[");
        String[] names={"minx","miny","minz","maxx","maxy","maxz"};
        int[][] arrays=new int[names.length][];
        for(int i=0;i<names.length;i++)arrays[i]=(int[])get(atlas,names[i]);
        for(int k=0;k<classes;k++){
            if(k>0)out.append(',');
            out.append('[').append(arrays[0][k]).append(',').append(arrays[1][k]).append(',').append(arrays[2][k]).append(',')
                .append(arrays[3][k]).append(',').append(arrays[4][k]).append(',').append(arrays[5][k]).append(']');
        }
        return out.append(']').toString();
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=5)throw new IllegalArgumentException("mem.bin shapes.bin nx ny nz required");
        int nx=Integer.parseInt(args[2]),ny=Integer.parseInt(args[3]),nz=Integer.parseInt(args[4]),classes=11;
        float[] flatMem=read(Path.of(args[0]),nx*ny*nz*classes);
        float[] flatShapes=read(Path.of(args[1]),nx*ny*nz*classes);
        float[][][][] mem=new float[nx][ny][nz][classes];
        float[][][][] shapes=new float[classes][nx][ny][nz];
        int q=0;
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++){
            mem[x][y][z][k]=flatMem[q];
            shapes[k][x][y][z]=flatShapes[q++];
        }
        MultipleToadDeformableAtlas atlas=new MultipleToadDeformableAtlas(classes);
        atlas.setImageInfo(nx,ny,nz,1.0f,1.0f,1.0f,0,0,0,0);
        set(atlas,"ntx",nx);set(atlas,"nty",ny);set(atlas,"ntz",nz);
        set(atlas,"rtx",1.0f);set(atlas,"rty",1.0f);set(atlas,"rtz",1.0f);
        set(atlas,"x0t",nx/2.0f);set(atlas,"y0t",ny/2.0f);set(atlas,"z0t",nz/2.0f);
        set(atlas,"rsx",new float[classes]);set(atlas,"rsy",new float[classes]);set(atlas,"rsz",new float[classes]);
        atlas.setShapes(shapes);
        for(int k=0;k<classes;k++)atlas.setShapeInfo(k,nx,ny,nz,1.0f,1.0f,1.0f);
        set(atlas,"registeredShape",new boolean[]{false,true,true,true,true,true,true,true,true,true,true});
        ParametricTransform model=new ParametricTransform("rigid",nx/2.0f,ny/2.0f,nz/2.0f,1,1,1,nx,ny,nz,
            nx/2.0f,ny/2.0f,nz/2.0f,1,1,1,nx,ny,nz);
        set(atlas,"transformModel",model);
        set(atlas,"transformMode",0); // NONE; updateShapeRegistration takes selected SINGLE route.
        atlas.setTransform(new float[]{0,0,0,0,0,0});
        byte[][][][] best=new byte[nx][ny][nz][4];
        byte[][][] segmentation=new byte[nx][ny][nz];
        atlas.updateShapeRegistration(mem,best,segmentation,new String[0],"rigid",5,1,0.0f,1.0f);
        atlas.registerShapes();
        float[] transform=atlas.getTransform();
        java.lang.reflect.Method energyMethod=atlas.getClass().getDeclaredMethod("computeExactRegistrationEnergy",float[].class);
        energyMethod.setAccessible(true);
        float finalEnergy=(Float)energyMethod.invoke(atlas,(Object)transform);
        System.out.println("{\"source\":\"VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + registerShapes() + private registerAllShapes()\","
            +"\"transform_bits\":"+vector(transform)+",\"energy_bits\":"+bits(finalEnergy)
            +",\"lambda_bits\":"+bits((Float)get(atlas,"lambda"))+",\"bounds\":"+bounds(atlas,classes)+"}");
    }
}
