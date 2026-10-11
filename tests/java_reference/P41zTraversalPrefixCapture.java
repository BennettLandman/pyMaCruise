import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Locale;

/** Capture one voxel after the in-place sweep prefix using the fixed Java oracle. */
public final class P41zTraversalPrefixCapture {
    private static Field field(Class<?> c, String name) throws Exception {
        Field f=c.getDeclaredField(name); f.setAccessible(true); return f;
    }
    private static Method method(Class<?> c, String name, Class<?>... args) throws Exception {
        Method m=c.getDeclaredMethod(name,args); m.setAccessible(true); return m;
    }
    private static String bits(float v) { return String.format(Locale.ROOT,"%08x",Float.floatToIntBits(v)); }
    private static String words(float[] a) {
        StringBuilder b=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append('"').append(bits(a[i])).append('"');}
        return b.append(']').toString();
    }
    public static void main(String[] args) throws Exception {
        Class<?> type=Class.forName("LesionNextIterationMembershipReferenceHarness");
        Object h=type.getDeclaredConstructor().newInstance();
        method(type,"read",String.class).invoke(h,args[0]);
        method(type,"computeSimilarities").invoke(h);
        int nx=(int)field(type,"nx").get(h), ny=(int)field(type,"ny").get(h), nz=(int)field(type,"nz").get(h);
        if(nx!=3 || ny!=30 || nz!=31) throw new IllegalArgumentException("expected cropped grid 3x30x31");
        int tx=1,ty=28,tz=29;
        Method update=method(type,"voxelUpdate",int.class,int.class,int.class);
        prefix: for(int x=1;x<nx-1;x++) for(int y=1;y<ny-1;y++) for(int z=1;z<nz-1;z++) {
            if(x==tx && y==ty && z==tz) break prefix;
            update.invoke(h,x,y,z);
        }
        int[] dx={1,-1,0,0,0,0},dy={0,0,1,-1,0,0},dz={0,0,0,0,1,-1};
        float[] memberships=(float[])field(type,"memberships").get(h);
        float[] lesion=(float[])field(type,"lesion").get(h);
        float[] black=(float[])field(type,"black").get(h);
        float[] image=(float[])field(type,"image").get(h);
        float[] shape=(float[])field(type,"shape").get(h);
        float[] mu=(float[])field(type,"mu").get(h), sd=(float[])field(type,"sd").get(h);
        float[] lesionMu=(float[])field(type,"lesionCentroid").get(h), lesionSd=(float[])field(type,"lesionStd").get(h);
        float[] weights=(float[])field(type,"weights").get(h);
        float[] shapePriors=(float[])field(type,"shapePriors").get(h), smoothingPriors=(float[])field(type,"smoothingPriors").get(h);
        int[] vd=(int[])field(type,"ventricleDistance").get(h), gd=(int[])field(type,"gmDistance").get(h), bd=(int[])field(type,"brainstemDistance").get(h);
        int n=nx*ny*nz, v=(tx*ny+ty)*nz+tz;
        float[] faceClasses=new float[6*13], faceLesion=new float[6], faceBlack=new float[6];
        for(int f=0;f<6;f++) {
            int q=((tx+dx[f])*ny+(ty+dy[f]))*nz+(tz+dz[f]);
            System.arraycopy(memberships,q*13,faceClasses,f*13,13);
            faceLesion[f]=lesion[q]; faceBlack[f]=black[q];
        }
        try(DataOutputStream out=new DataOutputStream(new FileOutputStream(args[1]))) {
            for(int c=0;c<2;c++) out.writeFloat(image[v*2+c]);
            for(float a:mu) out.writeFloat(a); for(float a:sd) out.writeFloat(a);
            for(float a:lesionMu) out.writeFloat(a); for(float a:lesionSd) out.writeFloat(a);
            for(float a:weights) out.writeFloat(a);
            for(int k=0;k<13;k++) out.writeFloat(shape[v*13+k]);
            for(float a:shapePriors) out.writeFloat(a); for(float a:smoothingPriors) out.writeFloat(a);
            for(float a:faceClasses) out.writeFloat(a); for(float a:faceLesion) out.writeFloat(a); for(float a:faceBlack) out.writeFloat(a);
            out.writeFloat(1.2f); out.writeFloat((0.2f/6.0f)/14.0f); out.writeFloat(2.0f/13.0f);
        }
        System.out.printf("{\"face_class_bits\":%s,\"face_lesion_bits\":%s,\"face_black_bits\":%s,\"target_shape10_bits\":\"%s\",",
            words(faceClasses),words(faceLesion),words(faceBlack),bits(shape[v*13+10]));
        update.invoke(h,tx,ty,tz);
        float[] lesionEnergy=(float[])field(type,"lesionEnergy").get(h);
        System.out.printf("\"target_final_lesion_energy_bits\":\"%s\"}%n",bits(lesionEnergy[v]));
    }
}
