import edu.jhmi.rad.medic.libraries.ImageFunctions;
import edu.jhmi.rad.medic.methods.ParametricTransform;
import java.io.DataInputStream;
import java.io.FileInputStream;

/** Source-selected ParametricTransform/ImageFunctions positive-support bounds. */
public final class LesionShapeBoundsReferenceHarness {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected input binary");
        try (DataInputStream in = new DataInputStream(new FileInputStream(args[0]))) {
            int nx=in.readInt(),ny=in.readInt(),nz=in.readInt();
            int tx=in.readInt(),ty=in.readInt(),tz=in.readInt();
            float[] ir={in.readFloat(),in.readFloat(),in.readFloat()};
            float[] tr={in.readFloat(),in.readFloat(),in.readFloat()};
            float[] transform=new float[6]; for(int i=0;i<6;i++)transform[i]=in.readFloat();
            float[][][] shape=new float[tx][ty][tz];
            for(int x=0;x<tx;x++)for(int y=0;y<ty;y++)for(int z=0;z<tz;z++)shape[x][y][z]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("unexpected trailing input");
            ParametricTransform model=new ParametricTransform("rigid",nx/2.0f,ny/2.0f,nz/2.0f,
                ir[0],ir[1],ir[2],nx,ny,nz,tx/2.0f,ty/2.0f,tz/2.0f,
                tr[0],tr[1],tr[2],tx,ty,tz);
            float[][] rotation=ParametricTransform.computeRotation(transform);
            int[] low={nx,ny,nz},high={0,0,0};int count=0;
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
                float[] point=model.imageToTemplate(x,y,z,transform,rotation,1.0f);
                if(ImageFunctions.linearInterpolation(shape,0.0f,point[0],point[1],point[2],tx,ty,tz)>0.0f) {
                    count++;if(x<low[0])low[0]=x;if(y<low[1])low[1]=y;if(z<low[2])low[2]=z;
                    if(x>high[0])high[0]=x;if(y>high[1])high[1]=y;if(z>high[2])high[2]=z;
                }
            }
            System.out.println("{\"lower_xyz\":["+low[0]+","+low[1]+","+low[2]+"],\"upper_xyz\":["+
                high[0]+","+high[1]+","+high[2]+"],\"support_count\":"+count+"}");
        }
    }
}
