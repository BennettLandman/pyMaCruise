import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import edu.jhmi.rad.medic.methods.MultipleToadDeformableAtlas;

/** Capture the selected single-scale rigid update at the Vandy iteration seam. */
public final class ToadsVandyPostRefreshRegistrationHarness {
    private static float[] read(Path path, int count) throws Exception {
        byte[] bytes=Files.readAllBytes(path);
        if(bytes.length!=count*4)throw new IllegalArgumentException("unexpected float payload: "+path);
        ByteBuffer b=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] out=new float[count];for(int i=0;i<count;i++)out[i]=b.getFloat();return out;
    }
    private static void set(Object target,String name,Object value)throws Exception{
        Field f=target.getClass().getDeclaredField(name);f.setAccessible(true);f.set(target,value);
    }
    private static Object get(Object target,String name)throws Exception{
        Field f=target.getClass().getDeclaredField(name);f.setAccessible(true);return f.get(target);
    }
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static String vector(float[] values){
        StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(bits(values[i]));}return s.append(']').toString();
    }
    private static String floats(float[] values){
        StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.toString(values[i]));}return s.append(']').toString();
    }
    public static void main(String[] args)throws Exception{
        if(args.length!=20)throw new IllegalArgumentException("expected membership/prior files, grid sizes, spacings, and six rigid parameters");
        int nx=Integer.parseInt(args[2]),ny=Integer.parseInt(args[3]),nz=Integer.parseInt(args[4]);
        int tx=Integer.parseInt(args[5]),ty=Integer.parseInt(args[6]),tz=Integer.parseInt(args[7]);
        float rix=Float.parseFloat(args[8]),riy=Float.parseFloat(args[9]),riz=Float.parseFloat(args[10]);
        float rtx=Float.parseFloat(args[11]),rty=Float.parseFloat(args[12]),rtz=Float.parseFloat(args[13]);
        float[] initial=new float[6];for(int i=0;i<6;i++)initial[i]=Float.parseFloat(args[14+i]);
        final int classes=11;
        float[] fm=read(Path.of(args[0]),nx*ny*nz*classes),fs=read(Path.of(args[1]),tx*ty*tz*classes);
        float[][][][] memberships=new float[nx][ny][nz][classes];
        float[][][][] shapes=new float[classes][tx][ty][tz];
        int q=0;
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)memberships[x][y][z][k]=fm[q++];
        q=0;
        for(int x=0;x<tx;x++)for(int y=0;y<ty;y++)for(int z=0;z<tz;z++)for(int k=0;k<classes;k++)shapes[k][x][y][z]=fs[q++];

        MultipleToadDeformableAtlas atlas=new MultipleToadDeformableAtlas(classes);
        set(atlas,"ntx",tx);set(atlas,"nty",ty);set(atlas,"ntz",tz);
        set(atlas,"rtx",rtx);set(atlas,"rty",rty);set(atlas,"rtz",rtz);
        set(atlas,"x0t",tx/2.0f);set(atlas,"y0t",ty/2.0f);set(atlas,"z0t",tz/2.0f);
        set(atlas,"nsx",new int[classes]);set(atlas,"nsy",new int[classes]);set(atlas,"nsz",new int[classes]);
        set(atlas,"rsx",new float[classes]);set(atlas,"rsy",new float[classes]);set(atlas,"rsz",new float[classes]);
        atlas.setImageInfo(nx,ny,nz,rix,riy,riz,0,0,0,0);
        atlas.setShapes(shapes);
        for(int k=0;k<classes;k++)atlas.setShapeInfo(k,tx,ty,tz,rtx,rty,rtz);
        boolean[] registered=new boolean[classes];for(int k=1;k<classes;k++)registered[k]=true;
        set(atlas,"registeredShape",registered);
        byte[][][][] best=new byte[nx][ny][nz][classes];
        byte[][][] segmentation=new byte[nx][ny][nz];
        atlas.initShapeRegistration(memberships,best,classes,"rigid",1,1);
        atlas.setTransform(initial.clone());
        atlas.updateShapeRegistration(memberships,best,segmentation,new String[0],"rigid",1,1,0.0f,1.0f);

        float[] hessian=new float[6],gradient=new float[6];
        Method coefficients=atlas.getClass().getDeclaredMethod("computeExactRegistrationCoefficients",float[].class,float[].class,float[].class);
        coefficients.setAccessible(true);
        float energy=(Float)coefficients.invoke(atlas,hessian,gradient,(Object)initial.clone());
        float[] sourceGradient=gradient.clone(),sourceHessian=hessian.clone();
        atlas.registerShapes();
        float[] finalTransform=atlas.getTransform().clone();
        Method energyMethod=atlas.getClass().getDeclaredMethod("computeExactRegistrationEnergy",float[].class);
        energyMethod.setAccessible(true);
        float finalEnergy=(Float)energyMethod.invoke(atlas,(Object)finalTransform);
        System.out.println("{\"initial_transform_bits\":"+vector(initial)+",\"energy_bits\":"+bits(energy)
            +",\"energy\":"+Float.toString(energy)+",\"gradient_bits\":"+vector(sourceGradient)
            +",\"gradient\":"+floats(sourceGradient)+",\"hessian_bits\":"+vector(sourceHessian)
            +",\"hessian\":"+floats(sourceHessian)+",\"final_transform_bits\":"+vector(finalTransform)
            +",\"final_energy_bits\":"+bits(finalEnergy)+",\"final_energy\":"+Float.toString(finalEnergy)
            +",\"lambda_bits\":"+bits((Float)get(atlas,"lambda"))+"}");
    }
}
