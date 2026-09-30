/*
 * Selected-method reference harness for VandyCRUISE TOADS dura-removal
 * memberships, scalar-image uncrop/import, hard-label byte packing, and the
 * selected MedicAlgorithmToads08 wrapper's unconditional ×255 scaling.
 *
 * The loops in generateDuraRemovalOutputs below are excerpted from
 * MultipleToadSegmentation.java; uncropAndBuffer is excerpted from
 * MultipleToadPreprocess.java. AlgorithmMultipleToads' selected output branch
 * calls that overload four times with a 0.0f background and imports each
 * buffer into one FLOAT image. The selected lower-case-option wrapper excerpt
 * has its scaleData parameter and guard commented out, then scales those four
 * ImageData outputs by 255.0f; the harness applies that selected operation to
 * the imported float buffers directly because the ImageDataMath dependency
 * source/runtime is absent from the available checkout. The wrapper also
 * receives the byte classification buffer with background 1. No MIPAV/JIST
 * image serialization is represented here. The CerebellarGM) misspelling is
 * intentionally kept.
 */
public final class ToadsDuraRemovalOutputsReferenceHarness {
    static final int NX=4, NY=3, NZ=3;
    static final int X0=-1, XN=1, Y0=0, YN=1, Z0=-1, ZN=0, BS=1;
    static final String[] NAMES={"SulcalCSF","Sulcal-CSF","CorticalGM",
        "Cerebellum-GM","CerebellarGM","CerebellarGM)","WhiteMatter","Background"};
    static final int K=NAMES.length;

    static float[][][][] memberships(int sx, int sy, int sz) {
        float[][][][] values=new float[sx][sy][sz][K];
        for (int x=0;x<sx;x++) for (int y=0;y<sy;y++) for (int z=0;z<sz;z++)
            for (int k=0;k<K;k++) values[x][y][z][k]=(float)(0.01*(1+x+10*y+100*z+1000*k));
        return values;
    }
    static int[][][] segmentation(int sx, int sy, int sz) {
        int[][][] values=new int[sx][sy][sz];
        for (int x=0;x<sx;x++) for (int y=0;y<sy;y++) for (int z=0;z<sz;z++)
            values[x][y][z]=(x+2*y+z)%K;
        return values;
    }

    // Source excerpt: MultipleToadSegmentation.generateDuraRemovalOutputs.
    static float[][][][] generate(float[][][][] mems, int[][][] segmentation) {
        int nx=segmentation.length, ny=segmentation[0].length, nz=segmentation[0][0].length;
        float[][][][] Mems=new float[4][nx][ny][nz];
        for (int k=0;k<K;k++) {
            if (NAMES[k].equals("SulcalCSF") || NAMES[k].equals("Sulcal-CSF"))
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
                    Mems[0][x][y][z]=mems[x][y][z][k];
        }
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
            Mems[1][x][y][z]=0.0f; Mems[2][x][y][z]=0.0f; Mems[3][x][y][z]=0.0f;
        }
        for (int k=0;k<K;k++) {
            if (!NAMES[k].equals("SulcalCSF") && !NAMES[k].equals("Sulcal-CSF")
                && !NAMES[k].equals("CorticalGM") && !NAMES[k].equals("Cerebrum-GM")
                && !NAMES[k].equals("CerebralGM") && !NAMES[k].equals("Cerebellum-GM")
                && !NAMES[k].equals("CerebellarGM") && !NAMES[k].equals("Background")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    Mems[2][x][y][z]+=(mems[x][y][z][k]);
                    if (segmentation[x][y][z]==k) Mems[3][x][y][z]=1.0f;
                }
            }
        }
        for (int k=0;k<K;k++) {
            if (NAMES[k].equals("CorticalGM") || NAMES[k].equals("Cerebrum-GM")
                || NAMES[k].equals("CerebralGM") || NAMES[k].equals("Cerebellum-GM")
                || NAMES[k].equals("CerebellarGM")) {
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    Mems[1][x][y][z]+=mems[x][y][z][k];
                    if (NAMES[segmentation[x][y][z]].equals("Cerebellum-GM")
                        || NAMES[segmentation[x][y][z]].equals("CerebellarGM)"))
                        Mems[3][x][y][z]=1.0f;
                }
            }
        }
        return Mems;
    }

    // Source excerpt: MultipleToadPreprocess.uncropAndBuffer(float[][][],float).
    static float[] uncrop(float[][][] src) {
        float[] larger=new float[NX*NY*NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            larger[x+NX*y+NX*NY*z]=0.0f;
        for (int x=Math.max(X0-BS,0);x<=Math.min(XN+BS,NX-1);x++)
            for (int y=Math.max(Y0-BS,0);y<=Math.min(YN+BS,NY-1);y++)
                for (int z=Math.max(Z0-BS,0);z<=Math.min(ZN+BS,NZ-1);z++)
                    larger[x+NX*y+NX*NY*z]=src[x-X0+BS][y-Y0+BS][z-Z0+BS];
        return larger;
    }
    // Source excerpt: MultipleToadPreprocess.uncropAndBuffer(byte[][][],byte).
    static byte[] uncrop(byte[][][] src) {
        byte[] larger=new byte[NX*NY*NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++)
            larger[x+NX*y+NX*NY*z]=(byte)1;
        for (int x=Math.max(X0-BS,0);x<=Math.min(XN+BS,NX-1);x++)
            for (int y=Math.max(Y0-BS,0);y<=Math.min(YN+BS,NY-1);y++)
                for (int z=Math.max(Z0-BS,0);z<=Math.min(ZN+BS,NZ-1);z++)
                    larger[x+NX*y+NX*NY*z]=src[x-X0+BS][y-Y0+BS][z-Z0+BS];
        return larger;
    }
    // Numeric operation selected by MedicAlgorithmToads08's ImageDataMath calls.
    static void scaleFloatValue(float[] values, float scale) {
        for (int i=0;i<values.length;i++) values[i]*=scale;
    }
    static void print(String name, float[] values, boolean comma) {
        if (comma) System.out.print(",");
        System.out.print("\""+name+"\":[");
        for (int i=0;i<values.length;i++) { if (i>0) System.out.print(","); System.out.print(Float.toString(values[i])); }
        System.out.print("]");
    }
    static void print(String name, byte[] values, boolean comma) {
        if (comma) System.out.print(",");
        System.out.print("\""+name+"\":[");
        for (int i=0;i<values.length;i++) { if (i>0) System.out.print(","); System.out.print(values[i]); }
        System.out.print("]");
    }
    public static void main(String[] args) {
        int cx=XN-X0+1+2*BS, cy=YN-Y0+1+2*BS, cz=ZN-Z0+1+2*BS;
        int[][][] labels=segmentation(cx,cy,cz);
        float[][][][] generated=generate(memberships(cx,cy,cz),labels);
        System.out.print("{\"original_shape\":["+NX+","+NY+","+NZ+"],\"crop_shape\":["+cx+","+cy+","+cz+"]");
        byte[][][] croppedClassification=new byte[cx][cy][cz];
        for (int x=0;x<cx;x++) for (int y=0;y<cy;y++) for (int z=0;z<cz;z++) croppedClassification[x][y][z]=(byte)labels[x][y][z];
        print("classification",uncrop(croppedClassification),true);
        String[] outputs={"csf","gray_matter","white_matter_fill","white_matter_mask"};
        for (int k=0;k<4;k++) {
            float[][][] volume=new float[cx][cy][cz];
            for (int x=0;x<cx;x++) for (int y=0;y<cy;y++) for (int z=0;z<cz;z++) volume[x][y][z]=generated[k][x][y][z];
            float[] importedBuffer=uncrop(volume);
            print("raw_"+outputs[k],importedBuffer,true);
            scaleFloatValue(importedBuffer,255.0f);
            print("scaled_"+outputs[k],importedBuffer,true);
        }
        System.out.println("}");
    }
}
